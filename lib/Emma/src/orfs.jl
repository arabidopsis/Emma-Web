#NCBI translation table 5, minus UUG start
const startcodon = biore"(ATT)|(ATC)|(ATA)|(ATG)|(GTG)|(TTG)"d
const stopcodon = biore"(TAG)|(TAA)|(AGA)|(AGG)"d


lengths = Dict{String, Int64}()
open("/home/pavel/Emma/CDS_median_lengths.txt", "r") do infile
    for line in readlines(infile)
        names = split(line, "\t")
        lengths[first(names)] = parse(Int64, last(names))
    end
end

#Dictionary constant of median lengths
const cdslengths = lengths

Codon = LongSubSeq{DNAAlphabet{4}}

function codonmatches(seq::CircularSequence, pattern)::Vector{Vector{Int32}}
    frames = [Int32[] for f in 1:3]
    for m in eachmatch(pattern, seq.sequence[1:seq.length+2])
        n_certain(matched(m)) < 3 && continue
        i::Int32 = m.captured[1]
        push!(frames[mod1(i, 3)], i)
    end
    return frames
end

function getcodons(seq::CircularSequence, pattern)
    positions = [Int32[] for f in 1:3]
    codons = Dict{Int32, Codon}()
    for m in eachmatch(pattern, seq.sequence[1:seq.length+2])
        n_certain(matched(m)) < 3 && continue
        i::Int32 = m.captured[1]
        push!(positions[mod1(i, 3)], i)
        codons[i] = matched(m)
    end
    return positions, codons
end

function getorfs!(writer::FASTA.Writer, id::AbstractString, genome::CircularSequence, strand::Char, starts::Vector{Vector{Int32}}, stops::Vector{Vector{Int32}}, minORF::Int)
    glength = length(genome)
    for (f, frame) in enumerate(starts)
        nextstop = 0
        for (s, start) in enumerate(frame)
            start < nextstop && continue
            nextstopidx = searchsortedfirst(stops[f], start+1)
            nextstop = first(stops[mod1(f - mod1(glength, 3), 3)]) #frame of next stop when wrapping depends on genome length
            if nextstopidx <= length(stops[f])
                nextstop = stops[f][nextstopidx]
            end
            circulardistance(start, nextstop, glength) < minORF && continue
            if nextstop < start; nextstop += glength; end
            translation = BioSequences.translate(genome.sequence[start:(nextstop-1)], code = ncbi_trans_table[2])
            translation[1] = AA_M
            write(writer, FASTA.Record(id * "*" * strand * "*" * string(start) * "-" * string(nextstop), translation))
        end
    end
end

function orfsearch(id::AbstractString, genome::CircularSequence, fstarts::Vector{Vector{Int32}}, fstops::Vector{Vector{Int32}},
    rstarts::Vector{Vector{Int32}}, rstops::Vector{Vector{Int32}}, minORF::Int)
    writer = open(FASTA.Writer, "tmp.orfs.fa")
    getorfs!(writer, id, genome, '+', fstarts, fstops, minORF)
    getorfs!(writer, id, reverse_complement(genome), '-', rstarts, rstops, minORF)
    close(writer)
    hmmpath = joinpath(emmamodels, "cds", "all_cds.hmm")
    cmd = `hmmsearch --domtblout tmp.domt $hmmpath tmp.orfs.fa`
    outfile = "tmp.hmmsearch.out"
    run(pipeline(cmd, stdout=outfile))
    return "tmp.domt"
end

struct HMMmatch
    orf::String
    strand::Char
    #a1::String
    #tlen::Int32
    query::String
    #a2::String
    qlen::Int
    evalue::Float64
    score::Float64
    #bias::Float64
    #num::Int
    #of::Int
    #cEvalue::Float64
    #iEvalue::Float64
    #dscore::Float64
    #dbias::Float64
    hmm_from::Int32
    hmm_to::Int32
    ali_from::Int32
    ali_to::Int32
    #env_from::Int
    #env_to::Int
    #acc::Float64
    #description::String
end

function parse_domt(file::String, glength::Integer)
    matches = FeatureMatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            startswith(line, "#") && continue
            #convert from ORF coordinates to genome coordinates (5'-3')
            bits = split(line, " ", keepempty=false)
            evalue = parse(Float64, bits[13]) # i-Evalue (indepent Evalue)
            evalue > 1 && continue # filter out poor matches
            orf = bits[1]
            orfbits = split(orf, "*")
            ends = split(orfbits[3], "-")
            strand = orfbits[2][1]
            orfstart = parse(Int32, ends[1])
            orfalifrom = parse(Int32, bits[18])
            orfalito = parse(Int32, bits[19])
            ali_from = orfstart + 3*(orfalifrom-1)
            ali_length = 3*(orfalito - orfalifrom + 1)
            if ali_from > glength
                ali_from = mod1(ali_from, glength)
            end
            
            # note that model coordinates are converted to nucleotide coordinates
            push!(matches, FeatureMatch(orf, bits[4], strand, 3 * parse(Int, bits[16]) - 2, 3 * parse(Int, bits[17]), ali_from, ali_length, evalue))
        end
    end
    @debug matches
    rationalise_matches!(matches, glength)
end

#const target_encoding = Dict(LongSequence{DNAAlphabet{2}}("ATT") => 0.0195429, LongSequence{DNAAlphabet{2}}("TTG") => 0.0280438, LongSequence{DNAAlphabet{2}}("ATG") => 0.180693,
#    LongSequence{DNAAlphabet{2}}("GTG") => 0.0286169, LongSequence{DNAAlphabet{2}}("ATA") => 0.024437, LongSequence{DNAAlphabet{2}}("ATC") => 0.0189099)

const target_encoding = Dict(LongSequence{DNAAlphabet{4}}("ATT") => 0.0047967167599361595, LongSequence{DNAAlphabet{4}}("TTG") => 0.0003434050417282311, LongSequence{DNAAlphabet{4}}("ATG") => 0.28243362577028663,
LongSequence{DNAAlphabet{4}}("GTG") => 0.04468146133661947, LongSequence{DNAAlphabet{4}}("ATA") => 0.018433691898543245, LongSequence{DNAAlphabet{4}}("ATC") => 0.0014133204427584738)

function fix_start_and_stop_codons!(hmm_matches, trns, starts, startcodons, stops, startcodon_model, glength)

    function is_possible_start(start, model_start, leftwindow, rightwindow, glength)
        d = closestdistance(start, model_start, glength) 
        #must be in frame
        mod(d, 3) ≠ 0 && return false
        if d > glength/2; d -= glength; end
        d > leftwindow && return false
        d < -rightwindow  && return false
        return true
    end

    sort!(hmm_matches; by=x->x.target_from)
    for (i,hmm_match) in enumerate(hmm_matches)
        @debug hmm_match
        hmmstart = hmm_match.target_from
        upstream_cds = hmm_matches[mod1(i-1, length(hmm_matches))]
        @debug "upstream_cds: $(upstream_cds.query)"
        #this works because if no trn gene is upstream of the CDS start, this will select the last trn gene in the genome, which is upstream of the first CDS
        trnidx = searchsortedfirst(trns, hmmstart, lt=(t,x)->t.fm.target_from < x)
        upstream_tRNA = trns[mod1(trnidx-1, length(trns))]
        @debug "upstream_tRNA: $(upstream_tRNA.fm.query)"
        inframe_stops = Int32[]
        for s in stops
            append!(inframe_stops, filter(x->mod(circulardistance(x, hmmstart, glength),3) == 0, s))
        end
        sort!(inframe_stops)
        stop_idx = searchsortedfirst(inframe_stops, hmmstart)
        upstream_stop = inframe_stops[mod1(stop_idx-1, length(inframe_stops))]
        distance_to_upstream_stop = circulardistance(upstream_stop, hmmstart, glength)
        @debug "$hmmstart, $distance_to_upstream_stop"
        #println(starts)
        possible_starts = Int32[]
        for s in starts
            append!(possible_starts, filter(x->is_possible_start(x, hmmstart, distance_to_upstream_stop, 50, glength), s))
        end
        if isempty(possible_starts)
            @warn "no starts for $hmm_match"
            continue
        end
        @debug possible_starts
        #Finds distance to upstream tRNA
        upstream_trn_end = upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length - 1
        dist_to_upstream_trn = closestdistance(upstream_trn_end, hmmstart, glength)
        @debug "dist to upstream tRNA: $dist_to_upstream_trn"
        #If distance is < 50, pick first inframe start of transcript
        if dist_to_upstream_trn < 10
            for (i,ps) in enumerate(possible_starts)
                relative_to_upstream_trn = closestdistance(upstream_trn_end, ps, glength)
                if relative_to_upstream_trn < 0
                    filter!(x -> x != ps, possible_starts)
                end
            end
            beststart = first(sort!(possible_starts))
            @debug "best start: $beststart"
        #If no nearby upstream tRNA, use XGboost
        else
            model_inputs = zeros(Float64, length(possible_starts), 6)
            for (i,ps) in enumerate(possible_starts)
                #calculate model inputs :target_encoding, :relative_to_hmm, :phase_to_hmm,:relative_to_upstream_tRNA, :relative_to_upstream_CDS, :relative_to_upstream_stop
                model_inputs[i, 1] = get(target_encoding, startcodons[ps], 0)
                relative_to_hmm = closestdistance(hmmstart, ps, glength)
                model_inputs[i, 2] = relative_to_hmm
                model_inputs[i, 3] = mod(relative_to_hmm, 3)
                upstream_trn_end = upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length - 1
                relative_to_upstream_trn = closestdistance(upstream_trn_end, ps, glength)
                model_inputs[i, 4] = relative_to_upstream_trn
                upstream_cds_end = upstream_cds.target_from + upstream_cds.target_length - 1
                relative_to_upstream_cds = closestdistance(upstream_cds_end, ps, glength)
                model_inputs[i, 5] = relative_to_upstream_cds
                relative_to_upstream_stop = closestdistance(upstream_stop, ps, glength)
                model_inputs[i, 6] = relative_to_upstream_stop
                @debug(model_inputs[i,:])
            end
            #predict with XGBoost model
            scores = XGBoost.predict(startcodon_model, model_inputs)
            @debug scores
            #pick top scoring start codon
            maxscore, maxidx = findmax(scores)
            beststart = possible_starts[maxidx]
            @debug "best start: $beststart $maxscore"
        end
        #pick first stop, or if none before tRNA, see if can construct stop by polyadenylation
        stop_idx = searchsortedfirst(stops[mod1(beststart, 3)], beststart) #index of first in-frame stop following best start codon
        next_stop = stops[mod1(beststart, 3)][mod1(stop_idx, length(stops[mod1(beststart, 3)]))]
        if closestdistance(beststart, next_stop, glength) > 2 * cdslengths[hmm_match.query]
            @warn "$(hmm_match.query) incomplete at the 3' boundary"
            next_stop = hmm_match.target_from + hmm_match.target_length - 1
        end
        @debug "first stop: $next_stop"
        #modify HMM match
        cds =  FeatureMatch(hmm_match.id, hmm_match.query, hmm_match.strand, 
                hmm_match.model_from, hmm_match.model_to, beststart, circulardistance(beststart, next_stop, glength), hmm_match.evalue)
        replace!(hmm_matches, hmm_match=>cds)
    end
    return hmm_matches
end
