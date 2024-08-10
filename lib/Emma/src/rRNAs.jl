const model2rrn = Dict("16srna" => "rrnL", "12srna" => "rrnS")

function parse_rrn_alignments(file::String, glength::Integer)
    alignments = FeatureMatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            !startswith(line, ">>") && continue
            readline(infile) #header line
            readline(infile) #dashes
            bits = split(readline(infile), " ", keepempty=false)
            readline(infile) #blank
            readline(infile) #NC
            readline(infile) #CS
            qseqline = strip(readline(infile))
            query = qseqline[1:(findfirst(" ", qseqline)[1]-1)]
            qseqline = lstrip(qseqline[length(query)+1:end])
            qseq = qseqline[(findfirst(" ", qseqline)[1]+1):(findlast(" ", qseqline)[1]-1)]
            readline(infile) #matches
            tseqline = strip(readline(infile))
            target = tseqline[1:(findfirst(" ", tseqline)[1]-1)]
            tstrand = bits[12][1]
            target_from = parse(Int, bits[10])
            tto = parse(Int, bits[11])
            if tstrand == '-'
                target_from = reverse_complement(target_from, glength)
                tto = reverse_complement(tto, glength)
            end
            rrn = model2rrn[query]
            qfrom = parse(Int, bits[7])
            push!(alignments, FeatureMatch(target, rrn, tstrand, qfrom, parse(Int, bits[8]), target_from, tto - target_from + 1, parse(Float64, bits[3])))
        end
    end
    return rationalise_matches!(alignments, glength)
end

function rrnsearch(tempfile::TempFile)
    hmmpath = joinpath(emmamodels, "rrn", "all_rrn.hmm")
    extended = tempfilename(tempfile, "tmp.extended.fa")
    tbl = tempfilename(tempfile, "tmp.tbl")
    cmd = `nhmmer --tblout $tbl $hmmpath $extended`

    outfile = tempfilename(tempfile, "tmp.nhmmer.out")
    run(pipeline(cmd, stdout=outfile))
    return tbl
end

function parse_tbl(file::String, glength::Integer)
    matches = FeatureMatch[]
    open(file, "r") do infile
        while !eof(infile)
            line = readline(infile)
            startswith(line, "#") && continue
            bits = split(line, " ", keepempty=false)
            rrnstrand = bits[12][1]
            rrnstart = parse(Int, bits[7])
            rrnstop = parse(Int, bits[8])
            fmlength = max(rrnstart, rrnstop) - min(rrnstart, rrnstop) + 1
            fmstart = rrnstrand == '+' ? mod1(rrnstart, glength) : mod1(reverse_complement(rrnstart, glength), glength)
            push!(matches, FeatureMatch(model2rrn[bits[3]], bits[3], rrnstrand, parse(Int, bits[5]), parse(Int, bits[6]),
                fmstart, fmlength, parse(Float64, bits[13])))
        end
    end
    @debug matches
    rationalise_matches!(matches, glength)
end

function find_closest_downstream_trn(target, trns, glength)
    idx = argmin(abs.([closestdistance(target, trn.fm.target_from + trn.fm.target_length, glength) for trn in trns]))
    closest_trn = trns[idx]
    return closest_trn
end

function fix_rrn_ends!(rRNAs, ftrns, rtrns, glength)
    for rrn in rRNAs
        trns = rrn.strand == '+' ? ftrns : rtrns
        rrnstart = rrn.target_from
        rrnstop = rrn.target_from + rrn.target_length - 1
        @debug "$rrnstart $rrnstop"
        changed = false
        trnidx = searchsortedfirst(trns, rrnstart, lt=(t, x) -> t.fm.target_from < x)
        upstream_tRNA = trns[mod1(trnidx - 1, length(trns))]
        @debug upstream_tRNA
        clockwise_dist = circulardistance(upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length, rrnstart, glength)
        anticlockwise_dist = circulardistance(rrnstart, upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length, glength)
        if clockwise_dist < 50 || anticlockwise_dist < rrnstart + rrn.target_length #overlap, or near enough to assume contiguity
            rrnstart = upstream_tRNA.fm.target_from + upstream_tRNA.fm.target_length
            changed = true
        end
        @debug "$rrnstart $rrnstop $changed"
        trnidx = searchsortedfirst(trns, rrnstop, lt=(t, x) -> t.fm.target_from + t.fm.target_length < x)
        #downstream_tRNA = trns[mod1(trnidx, length(trns))]
        downstream_tRNA = find_closest_downstream_trn(rrnstop, trns, glength)
        @debug downstream_tRNA
        if downstream_tRNA.fm.target_from - 50 < rrnstop #overlap, or near enough to assume contiguity
            rrnstop = downstream_tRNA.fm.target_from - 1
            changed = true
        end
        @debug "$rrnstart $rrnstop $changed"
        if changed
            replace!(rRNAs, rrn => FeatureMatch(rrn.id, rrn.query, rrn.strand, rrn.model_from, rrn.model_to, rrnstart,
                circulardistance(rrnstart, rrnstop, glength) + 1, rrn.evalue))
        end
    end
end
