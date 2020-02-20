# Segment cells with Julia

Welcome to SegmentationTools.jl! This is a random collection of functions that I've implemented for segmenting single cells from microscopy images.

```@meta
CurrentModule=SegmentationTools
```

## Flow

The general procedure follows the following format:

Segmentation -> Linking -> Analysis

## Segmentation

```@docs
SegmentationTools.segment_cells
```

## Linking

```@docs
SegmentationTools.build_tp_df
```

I generally use it like so:

```julia
particles = SegmentationTools.build_tp_df(img, thresholds);
using PyCall
tp = pyimport("trackpy")

import Pandas
# convert particles to a Pandas DataFrame and call trackpy
t = tp.link_df(Pandas.DataFrame(particles), 5);
linked = DataFrame(Pandas.DataFrame(t)); # convert back to Julian dataframe
```

## Analysis and Plotting

It's often helpful to look at each cell quickly to validate the tracking and
segmentation, which is where [`create_cell_grid`](@ref) comes in handy. Let's
create a XYT image that we've segmented and tracked.

The `tracks` dataframe should have the following columns (`trackpy` outputs
these by default):

* `x`, `y`: the centroid coordinates
* `frame`: the time index
* `id`: the per-frame id that can vary over the time course.
* `particle`: a persistent `id` that is the same across all time points for a
  given particle. `trackpy` links up the different `id`s into one `particle`
  value.

```@example cell_grid
using ImageShow
using DataFrames
using Colors
using SegmentationTools

img = zeros(Int, 25, 45, 1)
xs = 11:12:35
img[13, xs, 1] .= [1,2,3] # three equidistant dots
tracks = DataFrame(y=13, x=xs, frame=1:3, particle=1:3, id=1:3)
Gray.(img[:, :, 1]./4.0)
```


The corresponding cell grid will have the `particle` id displayed in the upper
right corner and show the time course of the movie with the particle's centroid
centered in the grid cell.

```@example cell_grid
result = SegmentationTools.create_cell_grid(img, tracks, win=10)
Gray.(result[:, :, 1])
```

```@docs
SegmentationTools.create_cell_grid
```

#### Observing cell localities

For low signal-to-noise applications, it can be useful to compute the local
background to normalize against. The exported function,
[`SegmentationTools.get_localities`](@ref), can be used to display the
localities used in [`SegmentationTools.build_tp_df`](@ref) to diagnose any weird
behaviors.


```@example
using FileIO, Images, SegmentationTools
slice = FileIO.load("https://user-images.githubusercontent.com/1661487/64392405-81a30500-d001-11e9-8262-744f76ff2f69.png")
labels = label_components(Bool.(slice))
foreground = labels .> 0.0

locals = Array{RGB{Float64}, 2}[]

gold = RGB((255,215,0)./255...)
orchid = RGB((186,85,211)./255...)

ids = sort(unique(labels))[2:end]
localities = SegmentationTools.get_localities(labels, ids, dist=(10, 60))

for id in ids
    mask = labels .== id
    result = RGB(0.5, 0.5, 0.5) .* foreground
    result[localities[id]] .= gold
    result[mask] .= orchid
    push!(locals, result)
end
reshape(locals, (:, 6))
```

```@docs
SegmentationTools.get_localities
```

## Miscellaneous

```@docs
SegmentationTools.find_peak_center
SegmentationTools.light_source_contrib
SegmentationTools.flatfield_correct
```