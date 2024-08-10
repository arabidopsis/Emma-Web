mutable struct GFF
    source::String
    ftype::String
    fstart::String
    fend::String
    score::String
    strand::Char
    phase::String
    attributes::String
end


function FMcoords2GFF(strand::Char, start::Integer, length::Integer, glength::Integer)
    gffstart = strand == '+' ? start : mod1(reverse_complement(start + length - 1, glength), glength)
    gffend = strand == '+' ? start + length - 1 : reverse_complement(start, glength)
    if gffend < gffstart
        gffend += glength
    end
    gffstart, gffend
end

function FMcoords2GFF(fm::FeatureMatch, glength::Integer)
    FMcoords2GFF(fm.strand, fm.target_from, fm.target_length, glength::Integer)
end

function CDS2GFF(cds::FeatureMatch, genome::CircularSequence, rev_genome::CircularSequence, trns::Vector{tRNA})
    glength = length(genome)
    name = first(split(cds.query, '.'))
    attributes = "Name=" * name * ";ID=$name.CDS"
    cdsstart = cds.target_from
    cdsstop = cds.target_from + cds.target_length - 1 + 3
    trnidx = findfirst(t -> circularin(t.fm.target_from, cdsstop - 2, 3, glength), trns)
    if !isnothing(trnidx)
        stopcodon = cds.strand == '+' ? copy(genome[cdsstop-2:cdsstop-2+circulardistance(cdsstop - 2, trns[trnidx].fm.target_from, glength)-1]) : copy(rev_genome[cdsstop-2:cdsstop-2+circulardistance(cdsstop - 2, trns[trnidx].fm.target_from, glength)-1])
        cdsstop = trns[trnidx].fm.target_from - 1
        while length(stopcodon) < 3
            push!(stopcodon, DNA_A)
        end
        attributes *= ";Note=putative $stopcodon stop codon is completed by the addition of 3' A residues to the mRNA"
    end
    gffstart, gffend = FMcoords2GFF(cds.strand, cdsstart, circulardistance(cdsstart, cdsstop + 1, glength), glength)
    return GFF("Emma", "CDS", string(gffstart), string(gffend), string(cds.evalue), cds.strand, "0", attributes)
end

function tRNA2GFF(trn::tRNA, glength::Integer; count="")
    gffstart, gffend = FMcoords2GFF(trn.fm, glength)
    attributes = "Name=" * trn.fm.query * "-" * trn.anticodon * ";ID=$(trn.fm.query)" * count
    if trn.polyA > 0
        attributes *= ";Note=tRNA completed by post-transcriptional addition of " * string(trn.polyA)
        attributes *= trn.polyA > 1 ? " As" : " A"
    end
    if trn.fm.query == "trnD" && trn.anticodon == "GCC"
        attributes *= ";Note=GUC anticodon completed by RNA editing"
        @warn "trnD $(trn.anticodon) edited to GUC by RNA editing"
    elseif !haskey(anticodon2trn, trn.anticodon)
        attributes *= ";Note=tRNA has no valid anticodon"
        @warn "no valid anticodon found for $(trn.fm.query)"
    end
    if typeof(attributes) != Missing
        return GFF("Emma", "tRNA", string(gffstart), string(gffend), string(trn.fm.evalue), trn.fm.strand, ".", attributes)
    else
        return nothing
    end
end

function rRNA2GFF(rrn::FeatureMatch, glength::Integer)
    gffstart, gffend = FMcoords2GFF(rrn, glength)
    attributes = "Name=" * rrn.query * ";ID=$(rrn.query)"
    return GFF("Emma", "rRNA", string(gffstart), string(gffend), string(rrn.evalue), rrn.strand, ".", attributes)
end

function add_geneGFFs(gffs::Vector{GFF}, genome::CircularSequence, glength::Integer)
    cdsGFFs = filter(x -> x.ftype == "CDS", gffs)
    for gff in cdsGFFs
        name = match(r"Name=(\w+)", gff.attributes).captures[1]
        #Checks for frameshift, to split feature into 2 CDS (note: code repeated from features.jl)
        if name == "ND3"
            notes = split(gff.attributes, ';')
            length(notes) < 3 ? notes = "" : notes = join(notes[3:end])
            fstart = parse(Int32, gff.fstart)
            fend = parse(Int32, gff.fend)
            gene_seq = string(genome[fstart:fend])
            sequence_match = match(r".TT.CT.AGTAGC", gene_seq)
            if sequence_match !== nothing
                CDS1end = fstart + sequence_match.offset + 4
                CDS2start = fstart + sequence_match.offset + 6
                deleteat!(gffs, findfirst(x -> x == gff, gffs))
                push!(gffs, GFF("Emma", "CDS", gff.fstart, string(CDS1end), gff.score, gff.strand, ".", "Name=$name;ID=$name.CDS.1;Note=CDS contains a frameshift where one nucleotide is skipped" * notes))
                push!(gffs, GFF("Emma", "CDS", string(CDS2start), gff.fend, gff.score, gff.strand, ".", "Name=$name;ID=$name.CDS.2;Note=CDS contains a frameshift where one nucleotide is skipped" * notes))
            end
        end
        #Adds gene features, that are essentially identical to CDS
        push!(gffs, GFF("Emma", "gene", gff.fstart, gff.fend, gff.score, gff.strand, ".", "Name=$name;ID=$name.gene"))
    end
    return gffs
end

function getGFF(genome::CircularSequence, rev_genome::CircularSequence, cds_matches::Vector{FeatureMatch},
    trn_matches::Vector{tRNA}, rRNAs::Vector{FeatureMatch}, glength::Integer)

    gffs = GFF[]
    genome_length = length(genome)

    grouped_trns = Dict{String,Vector{tRNA}}()
    for trn in trn_matches
        if haskey(grouped_trns, trn.fm.query)
            push!(grouped_trns[trn.fm.query], trn)
        else
            grouped_trns[trn.fm.query] = [trn]
        end
    end

    function writeone(gff::Union{Nothing,GFF})
        if typeof(gff) != Nothing
            push!(gffs, gff)
        end
    end

    for cds in cds_matches
        gff = CDS2GFF(cds, genome, rev_genome, trn_matches)
        writeone(gff)
    end
    for (id, trns) in grouped_trns
        if length(trns) > 1
            for (i, trn) in enumerate(trns)
                gff = tRNA2GFF(trn, genome_length, count=".$i")
                writeone(gff)
            end
        else
            gff = tRNA2GFF(trns[1], genome_length)
            writeone(gff)
        end
    end
    for rrn in rRNAs
        gff = rRNA2GFF(rrn, genome_length)
        writeone(gff)
    end
    #gffs = trnF_start(gffs, glength)
    gffs = add_geneGFFs(gffs, genome, glength)
    fcds = filter(x -> x.strand == '+', cds_matches)
    ftrns = filter(x -> x.fm.strand == '+', trn_matches)
    fmRNAs = add_mRNAGFFs(fcds, ftrns, glength)
    rcds = filter(x -> x.strand == '-', cds_matches)
    rtrns = filter(x -> x.fm.strand == '-', trn_matches)
    rmRNAs = add_mRNAGFFs(rcds, rtrns, glength)
    append!(gffs, fmRNAs, rmRNAs)
    return gffs
end

function writeGFF(id::AbstractString, gffs::Vector{GFF}, outfile::String)
    open(outfile, "w") do out
        writeGFF(id, gffs, out)
    end
end

function writeGFF(id::AbstractString, gffs::Vector{GFF}, outfile::IO)
    for gff in gffs
        write(outfile, join([id, gff.source, gff.ftype, gff.fstart, gff.fend, gff.score, gff.strand, gff.phase, gff.attributes], "\t"), "\n")
    end

end



function add_mRNAGFFs(cds_matches::Vector{FeatureMatch}, trns::Vector{tRNA}, glength::Integer)
    mRNAs = GFF[]
    for (i, cds_match) in enumerate(cds_matches)
        name = first(split(cds_match.query, '.'))
        attributes = "Name=" * name * ";ID=$name.mRNA"
        cdsstart = cds_match.target_from
        cdsstop = cdsstart + cds_match.target_length - 1 + 3
        trnidx = searchsortedfirst(trns, cdsstart, lt=(t, x) -> t.fm.target_from < x)

        upstream_tRNA = trns[mod1(trnidx - 1, length(trns))]
        upstream_trn_end = upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length - 1
        dist_to_upstream_trn = abs(closestdistance(upstream_trn_end, cdsstart, glength))

        downstream_tRNA = trns[mod1(trnidx, length(trns))]
        downstream_trn_start = downstream_tRNA.fm.target_from
        dist_to_downstream_trn = abs(closestdistance(downstream_trn_start, cdsstop, glength))
        if dist_to_upstream_trn < 10
            cdsstart = upstream_trn_end + 1
        else
            attributes *= ";Note=5' incomplete"
        end
        if dist_to_downstream_trn < 10
            cdsstop = downstream_trn_start - 1
        else
            attributes *= ";Note=3' incomplete"
        end
        gffstart, gffend = FMcoords2GFF(cds_match.strand, cdsstart, circulardistance(cdsstart, cdsstop + 1, glength), glength)
        push!(mRNAs, GFF("Emma", "mRNA", string(gffstart), string(gffend), string(cds_match.evalue), cds_match.strand, "0", attributes))
    end
    return mRNAs
end

