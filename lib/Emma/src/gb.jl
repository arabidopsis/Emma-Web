
products = Dict{String,String}()

open(joinpath(DATA, "gene2product.txt"), "r") do infile
    for line in readlines(infile)
        names = split(line, "\t")
        products[first(names)] = last(names)
    end
end

gene2products = products


function feature_compare(x::GFF, y::GFF)
    feature_hierarchy = Dict("gene" => 1, "CDS" => 2, "mRNA" => 3, "tRNA" => 4, "rRNA" => 5)
    return feature_hierarchy[x.ftype] < feature_hierarchy[y.ftype]
end

function writeGB(id::AbstractString, gffs::Vector{GFF}, outfile_gb::String)
    open(outfile_gb, "w") do out
        writeGB(id, gffs, out)
    end
end

function writeGB(id::AbstractString, gffs::Vector{GFF}, out::IO)
    write(out, ">Feature $id\n")
    sort!(gffs, by=x -> parse.(Int, x.fstart))
    gffs = group_features(gffs) #Ensures that features of common locus are written together
    for (name, vector) in gffs
        sort!(vector, lt=feature_compare) #Ensures correct heirarchy in file
        for gff in vector
            gff.strand == '+' ? fstart = gff.fstart : fstart = gff.fend
            gff.strand == '+' ? fend = gff.fend : fend = gff.fstart
            attributes = eachmatch(r"Note=([^;]+)", gff.attributes) #Matches notes and adds them to attributes vector
            if gff.ftype == "CDS"
                if "Note=CDS contains a frameshift where one nucleotide is skipped" in [m.match for m in attributes] #For special frameshift cases
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                    cdsend = parse(Int64, fend)
                    second_cds = vector[findfirst(x -> x.fstart == string(cdsend + 2), vector)]
                    write(out, join([second_cds.fstart, second_cds.fend], '\t'), '\n')
                    write(out, '\t'^3, "exception\tribosomal slippage\n")
                    deleteat!(vector, findfirst(x -> x == second_cds, vector))
                else
                    write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                end
                write(out, '\t'^3, "gene\t$name", '\n')
                write(out, '\t'^3, "product\t$(gene2products[name])", '\n')
                write(out, '\t'^3, "transl_table\t2", '\n')
                if !isempty(attributes)
                    note = join([last(split(m.match, '=')) for m in attributes], ", ") #Appends notes into single string
                    write(out, '\t'^3, "note\t$note", '\n')
                end

            elseif gff.ftype == "tRNA"
                anticodon = match(r"trn..-(\w+)", gff.attributes) #Matches trnS/L features to add anticodon note
                write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                write(out, '\t'^3, "product\t$(gene2products[name])", '\n')
                if anticodon != nothing
                    write(out, '\t'^3, "note\tAnticodon: $(anticodon.captures[1])", '\n')
                end

            elseif gff.ftype == "rRNA"
                write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                write(out, '\t'^3, "product\t$(gene2products[name])", '\n')

            elseif gff.ftype == "gene"
                write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                write(out, '\t'^3, "gene\t$name", '\n')

            elseif gff.ftype == "mRNA"
                note = nothing
                if !isempty(attributes)
                    if "Note=5' incomplete" in [m.match for m in attributes]
                        fstart = '<' * fstart
                    end
                    if "Note=3' incomplete" in [m.match for m in attributes]
                        fend = '>' * fend
                    end
                    note = join([last(split(m.match, '=')) for m in attributes], ", ")
                end
                write(out, join([fstart, fend, gff.ftype], '\t'), '\n')
                write(out, '\t'^3, "gene\t$name", '\n')
                write(out, '\t'^3, "product\t$(gene2products[name])", '\n')
                if note !== nothing
                    write(out, '\t'^3, "note\t$note", '\n')
                end
            end
        end
    end
end

#function write

function group_features(gffs::Vector{GFF})
    gff_dict = Dict{String,Vector{GFF}}()
    for gff in gffs
        name = match(r"Name=(\w+)", gff.attributes).captures[1]
        if haskey(gff_dict, name)
            push!(gff_dict[name], gff)
        else
            gff_dict[name] = [gff]
        end
    end
    return gff_dict
end
