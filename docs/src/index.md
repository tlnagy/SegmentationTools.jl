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
particles = SegmentationTools.build_tp_df(img1, segments, :EPI_mNG);
using PyCall
tp = pyimport("trackpy")

import Pandas
# convert particles to a Pandas DataFrame and call trackpy
t = tp.link_df(Pandas.DataFrame(particles), 5);
linked = DataFrame(Pandas.DataFrame(t)); # convert back to Julian dataframe
```

## Analysis and Plotting

```@example
using ImageShow
using DataFrames
using Colors
using SegmentationTools

img = zeros(Int, 25, 45, 1)
xs = 11:12:35
img[13, xs, 1] .= [1,2,3] # three equidistant dots
tracks = DataFrame(y=13, x=xs, frame=1:3, particle=1:3, id=1:3)
result = SegmentationTools.create_cell_grid(img, tracks, win=10)
Gray.(result[:, :, 1])
```

```@docs
SegmentationTools.create_cell_grid
```

## Miscellaneous

```@docs
SegmentationTools.find_peak_center
SegmentationTools.light_source_contrib
SegmentationTools.flatfield_correct
```