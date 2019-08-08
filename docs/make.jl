# SegmentationTools isn't published yet so we need to add the source directly
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Documenter, SegmentationTools

makedocs(sitename="SegmentationTools.jl")

