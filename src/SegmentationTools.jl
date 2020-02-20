module SegmentationTools

using AxisArrays
using ColorSchemes
using Colors
using DataFrames
using DataStructures
using FreeTypeAbstraction
using ImageAxes
using ImageSegmentation
using ImageMorphology
using ImageMetadata
using Images
using KernelDensity
using MappedArrays
using Statistics
using SimpleTraits

include("utils.jl")
include("segment.jl")
include("diagnostics.jl")

export get_localities

end # module
