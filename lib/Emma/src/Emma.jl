
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
using Pkg.Artifacts

export main, doone, writeGFF, tempfilename, TempFile, drawgenome, trnF_start


const DATA = joinpath(@__DIR__, "..")
const emmamodels = joinpath(DATA, "emma_vertebrate_models")
# const emmamodels = joinpath(artifact"Emma-models","emma-models-0.1.5-alpha")
# const emmamodels = joinpath(artifact"Emma2-models","emma-models-1.0.0")
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
