
module Emma
using Serialization
using Artifacts
using BioSequences
using FASTX
using XGBoost
using Logging
using Unicode
using GenomicAnnotations
using UUIDs

export main, doone, writeGFF, tempfilename, TempFile, drawgenome, trnF_start


const DATA = joinpath(@__DIR__, "..")
const emmamodels = joinpath(DATA, "emma_vertebrate_models")

include("tempfile.jl")
include("circularity.jl")
include("feature.jl")
include("orfs.jl")
include("tRNAs.jl")
include("rRNAs.jl")
include("gff.jl")
include("gb.jl")
include("visuals.jl")
include("process.jl")
include("cmd.jl")

end
