module Chloe2

#using Serialization
using ArgMacros
using Artifacts
using BioSequences
using FASTX
using Logging
using Unicode
using GenomicAnnotations
using UUIDs
using CSV
using DataFrames
using Printf

export main, chloe, chloeone, writeGFF, tempfilename, TempFile

const chloe2models = joinpath(artifact"Chloe2_models", "Chloe2_models-2.0.0-alpha.2", "models")

include("tempfile.jl")
include("circularity.jl")
include("feature.jl")
include("orfs.jl")
include("introns.jl")
include("tRNAs.jl")
include("rRNAs.jl")
include("annotations.jl")
include("gff.jl")
#include("gb.jl")
include("process.jl")
include("cmd.jl")

end
