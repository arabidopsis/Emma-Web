
module Emma
using Serialization
using Artifacts
using BioSequences
using FASTX
using XGBoost
using Logging
using Unicode
using GenomicAnnotations
using ArgMacros
using UUIDs

export main, doone, writeGFF


const DATA = joinpath(@__DIR__, "..")
const emmamodels = joinpath(DATA, "emma_vertebrate_models")
include("circularity.jl")
include("feature.jl")
include("orfs.jl")
include("tRNAs.jl")
include("rRNAs.jl")
include("gff.jl")
include("gb.jl")
include("visuals.jl")

const minORF = 150


function remove_starts_in_tRNAs!(starts, codons, tRNAs, glength)
    fms = getproperty.(tRNAs, :fm)
    for (startvector, codonvector) in zip(starts, codons)
        startsintRNA = []
        for (i, s) in enumerate(startvector)
            if any(circularin.(s, fms, glength))
                push!(startsintRNA, i)
            end
        end
        deleteat!(startvector, startsintRNA)
    end
end

function add_stops_at_tRNAs!(stops, tRNAs, glength)
    for t in tRNAs
        push!(stops[mod1(t.fm.target_from, 3)], t.fm.target_from)
        minus1 = mod1(t.fm.target_from - 1, glength)
        push!(stops[mod1(minus1, 3)], minus1)
        minus2 = mod1(t.fm.target_from - 2, glength)
        push!(stops[mod1(minus2, 3)], minus2)
    end
    sort!(stops[1])
    sort!(stops[2])
    sort!(stops[3])
    return stops
end

#If duplicate features exist, remove the ones with higher e-values
function remove_duplicate_features(matches::Vector)
    features = FeatureMatch[]
    sorted_array = sort(matches, by=x -> x.evalue)
    filtered_dict = Dict()
    for item in sorted_array
        if !haskey(filtered_dict, item.query)
            filtered_dict[item.query] = item
        end
    end
    filtered_array = collect(values(filtered_dict))
    for feature in filtered_array
        push!(features, feature)
    end
    return features
end

function trnF_start(GFFs, genome, glength)
    trnF_idx = findfirst(x -> occursin(r"trnF", x.attributes), GFFs)
    if trnF_idx != nothing
        trnF = GFFs[trnF_idx]
        offset = parse(Int32, trnF.fstart) - 1
        for gff in GFFs
            fstart = parse(Int32, gff.fstart)
            fend = parse(Int32, gff.fend)
            fstart -= offset
            if fstart < 1
                fstart += glength
            end
            fend -= offset
            if fend < 1
                fend += glength
            end
            gff.fstart = string(fstart)
            gff.fend = string(fend)
        end
        if offset > 1
            genome = LongDNA{4}(string(genome[offset+1:glength]) * string(genome[1:offset]))
        else
            genome = LongDNA{4}(string(genome[1:glength]))
        end
        return GFFs, genome
    else
        @warn "Positional translation not possible due to missing trnF"
        return GFFs, genome
    end
end

function cleanfiles(uuid::UUID)
    for f in ["tmp.cmsearch.out", "tmp.extended.fa", "tmp.nhmmer.out", "tmp.tbl",
        "tmp.domt", "tmp.hmmsearch.out", "tmp.orfs.fa"]
        rm("$(uuid).$(f)", force=true)
    end
end



function doone(infile::String, uuid::UUID)
    @info "$infile"


    target = FASTA.Record()
    reader = open(FASTA.Reader, infile)
    read!(reader, target)
    id = FASTA.identifier(target)
    genome = CircularSequence(FASTA.sequence(LongDNA{4}, target))
    glength = length(genome)
    rev_genome = reverse_complement(genome)
    @info "$id\t$glength bp"

    #extend genome
    extended_genome = genome[1:glength+100]
    writer = open(FASTA.Writer, "$(uuid).tmp.extended.fa")
    write(writer, FASTA.Record(id, extended_genome))
    close(writer)

    #find tRNAs
    trn_matches = parse_trn_alignments(cmsearch(uuid, "trn", "all_trn.cm"), glength)
    @debug trn_matches
    filter!(x -> x.fm.evalue < 1e-5, trn_matches)
    filter!(x -> x.fm.target_from <= glength, trn_matches)
    ftrns = get_best_trns(filter(x -> x.fm.strand == '+', trn_matches), glength)
    rtrns = get_best_trns(filter(x -> x.fm.strand == '-', trn_matches), glength)
    #check for overlapping tRNAs
    overlapped = get_overlapped_trns(sort(ftrns, by=x -> x.fm.target_from), glength)
    append!(overlapped, get_overlapped_trns(sort(rtrns, by=x -> x.fm.target_from), glength))
    trn_matches = append!(ftrns, rtrns)
    #for overlapped tRNAs, generate polyadenylated version
    for (cma, trunc_end) in overlapped
        trnseq = cma.fm.strand == '+' ? genome.sequence[cma.fm.target_from:trunc_end] : rev_genome.sequence[cma.fm.target_from:trunc_end]
        trnseq_polyA = trnseq * dna"AAAAAAAAAA"
        polyA_matches = parse_trn_alignments(cmsearch(uuid, cma.fm.query, "trn", trnseq_polyA), 0)
        isempty(polyA_matches) && continue
        trn_match = polyA_matches[1]
        if trn_match.fm.evalue < cma.fm.evalue
            #construct modified CMA with new evalue, new end point, polyA length
            newcma = tRNA(FeatureMatch(cma.fm.id, cma.fm.query, cma.fm.strand, trn_match.fm.model_from, trn_match.fm.model_to, cma.fm.target_from,
                    circulardistance(cma.fm.target_from, trunc_end, glength) + 1, trn_match.fm.evalue),
                trn_match.anticodon, (trn_match.fm.target_length) - (trunc_end - cma.fm.target_from))
            #delete old match from trn_matches
            deleteat!(trn_matches, findfirst(x -> x == cma, trn_matches))
            #add new CMA
            push!(trn_matches, newcma)
        end
    end

    @debug trn_matches
    @info "found $(length(trn_matches)) tRNA genes"

    #find rRNAs
    #search for rrns using hmmsearch
    rrns = parse_tbl(rrnsearch(uuid), glength)
    @debug rrns
    #fix ends using flanking trn genes
    ftRNAs = filter(x -> x.fm.strand == '+', trn_matches)
    sort!(ftRNAs, by=x -> x.fm.target_from)
    rtRNAs = filter(x -> x.fm.strand == '-', trn_matches)
    sort!(rtRNAs, by=x -> x.fm.target_from)
    fix_rrn_ends!(rrns, ftRNAs, rtRNAs, glength)
    @debug rrns

    #find CDSs
    fstarts, fstartcodons = getcodons(genome, startcodon)
    remove_starts_in_tRNAs!(fstarts, fstartcodons, ftRNAs, glength)
    fstops = codonmatches(genome, stopcodon)
    add_stops_at_tRNAs!(fstops, ftRNAs, glength)
    @debug fstops

    rstarts, rstartcodons = getcodons(rev_genome, startcodon)
    remove_starts_in_tRNAs!(rstarts, rstartcodons, rtRNAs, glength)
    rstops = codonmatches(rev_genome, stopcodon)
    add_stops_at_tRNAs!(rstops, rtRNAs, glength)

    cds_matches = parse_domt(orfsearch(uuid, id, genome, fstarts, fstops, rstarts, rstops, minORF), glength)
    @debug cds_matches
    #fix start & stop codons
    #load XGBoost model
    startcodon_model = Booster(DMatrix[], model_file=joinpath(emmamodels, "xgb.model"))
    fhmms = filter(x -> x.strand == '+', cds_matches)
    rhmms = filter(x -> x.strand == '-', cds_matches)
    fix_start_and_stop_codons!(fhmms, ftRNAs, fstarts, fstartcodons, fstops, startcodon_model, glength)
    fix_start_and_stop_codons!(rhmms, rtRNAs, rstarts, rstartcodons, rstops, startcodon_model, glength)

    cds_matches = append!(fhmms, rhmms)
    frameshift_merge!(cds_matches, glength, genome)
    cds_matches = remove_duplicate_features(cds_matches)
    @info "found $(length(cds_matches)) protein-coding genes"

    gffs = getGFF(genome, rev_genome, cds_matches, trn_matches, rrns, glength)
    cleanfiles(uuid)
    return id, gffs, genome
end
function main(infile::String; outfile_gff="", outfile_gb="", outfile_fa="", svgfile="", loglevel="Info")
    global_logger(ConsoleLogger(loglevel == "debug" ? Logging.Debug : Logging.Info))

    uuid = uuid4()
    id, gffs, genome = doone(infile, uuid)

    glength = length(genome)
    CDSless = filter(x -> x.ftype != "CDS", gffs)
    mRNAless = filter(x -> x.ftype != "mRNA", CDSless)
    if outfile_fa !== nothing
        gffs, shifted_genome = trnF_start(gffs, genome, glength)
        open(FASTA.Writer, outfile_fa) do w
            write(w, FASTA.Record(id, shifted_genome))
        end
    end
    if outfile_gff !== nothing
        writeGFF(id, gffs, outfile_gff)
    else
        for gff in gffs
            println(join([id, gff.source, gff.ftype, gff.fstart, gff.fend, gff.score, gff.strand, gff.phase, gff.attributes], "\t"))
        end
    end
    if outfile_gb !== nothing
        writeGB(id, gffs, outfile_gb, glength)
    end
    if svgfile !== nothing
        drawgenome(svgfile, id, glength, mRNAless)
    end
    # cleanfiles(uuid)
end



function main()
    args = @dictarguments begin
        @helpusage "Emma.jl [options] <FASTA_file>"
        @helpdescription """
            Note: Use consistant inputs/outputs. If you wish
            to annotate a directory of fasta files, ensure that
            the output parameters are also directories.
            """
        @argumentoptional String GFF_out "--gff"
        @arghelp "file/dir for gff output"
        @argumentoptional String GB_out "--gb"
        @arghelp "file/dir for gb output"
        @argumentoptional String FA_out "--fa"
        @arghelp "file/dir for fasta output. Use this argument if you wish annotations to begin with tRNA-Phe"
        @argumentoptional String Svg_out "--svg"
        @arghelp "file/dir for SVG output"
        @positionalrequired String FASTA_file
        @arghelp "file/dir for fasta input"
    end
    filtered_args = filter(pair -> pair.second !== nothing, args)
    all_dirs = all(isdir, values(filtered_args))

    function maybefile(t)
        isdir(t) && return false
        return true
    end
    all_files = all(maybefile, values(filtered_args))
    if all_dirs
        fafiles = filter!(x -> endswith(x, ".fa") || endswith(x, ".fasta"), readdir(args[:FASTA_file], join=true))
        for fasta in fafiles
            accession = first(split(basename(fasta), "."))
            outfile_gff = haskey(filtered_args, :GFF_out) ? joinpath(args[:GFF_out], accession * ".gff") : nothing
            outfile_gb = haskey(filtered_args, :GB_out) ? joinpath(args[:GB_out], accession * ".gb") : nothing
            outfile_fa = haskey(filtered_args, :FA_out) ? joinpath(args[:FA_out], accession * ".fa") : nothing
            outfile_svg = haskey(filtered_args, :Svg_out) ? joinpath(args[:Svg_out], accession * ".svg") : nothing
            main(fasta, outfile_gff=outfile_gff, outfile_gb=outfile_gb, outfile_fa=outfile_fa, svgfile=outfile_svg)
        end
    elseif all_files
        main(args[:FASTA_file], outfile_gff=args[:GFF_out], outfile_gb=args[:GB_out], outfile_fa=args[:FA_out], svgfile=args[:Svg_out])
    else
        throw("Inputs must be consistant; all directories or all files")
    end
end

end

#ARGS[1] = fasta input
#ARGS[2] = gff output
#ARGS[3] = svg output
