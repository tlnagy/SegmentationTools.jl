using Images: otsu_threshold, imadjustintensity
using Unitful
using Unitful: μm
using OffsetArrays
using ProgressMeter

"""
    segment_cells(img, seed_channel, segment_channel, maskfunc)
"""
function segment_cells(img::AxisArray{T, 4},
                       seed_channel::Symbol,
                       segment_channel::Symbol,
                       maskfunc::Function) where {T}
    if !all(map(x->x in axisnames(img), [:time, :channel]))
        error("image requires :channel and :time axes")
    end
    tax = timeaxis(img)
    segments = Array{SegmentedImage{Array{Int64,2},Float32}}(undef, length(tax))
    ylen, xlen = size(img, Axis{:y}), size(img, Axis{:x})
    mask = Array{Bool}(undef, ylen, xlen)
    segslice = Array{T}(undef, ylen, xlen)
    seedslice = Array{T}(undef, ylen, xlen)
    seeds = Array{Bool}(undef, ylen, xlen)

    for (idx, tp) in enumerate(tax.val)

        # find seeds
        seedslice .= view(img, Axis{:channel}(seed_channel), Axis{:time}(tp))
        seeds .= seedslice .> otsu_threshold(seedslice)*1.1

        # find segmentation slice
        segslice .= view(img, Axis{:channel}(segment_channel), Axis{:time}(tp))

        # throw out seeds that are outside of segment areas
        seeds .&= (segslice .> 0.0)
        ImageMorphology.closing!(seeds)
        ImageMorphology.erode!(seeds)
        markers = label_components(seeds)

        # compute mask
        mask .= maskfunc(segslice)
        ImageMorphology.closing!(mask)

        segments[idx] = watershed(1 .- imadjustintensity(segslice), markers, mask=mask)
    end
    segments
end

"""
    build_tp_df(img, components; dist)

Build a `DataFrames.DataFrame` that is compatible with trackpys `link_df`
function. Needs to be converted to a `Pandas.DataFrame` before passing to
trackpy. `dist` is a 2-tuple of integers indicating the minimum and maximum
distance away in pixels from each cell to include in its local background
calculation.
"""
function build_tp_df(img::AxisArray{T1, 4},
                     components::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Integer}

    xstep = step(AxisArrays.axes(img, Axis{:x}).val)
    ystep = step(AxisArrays.axes(img, Axis{:y}).val)
    (xstep != ystep) && @warn "Different scaling for x and y axes is not supported"
    pixelarea = xstep * ystep
    particles = DataFrames.DataFrame[]

    @showprogress 1 "Computing..." for (idx, timepoint) in enumerate(timeaxis(img))
        component_slice = view(components, :, :, idx)

        # get all ids present in the current slice
        slice_ids = sort(unique(component_slice))
        # get the number of pixels for each id in the slice
        lengths = Images.component_lengths(component_slice)[slice_ids .+ 1]

        if eltype(pixelarea) <: Unitful.Area
            # convert boolean area to microns with the assumption that the pixel
            # space is the same in both x and y
            areas = round.(μm^2, lengths .* pixelarea, sigdigits=4)
            # filter out too large and too small particles
            correct_size = (10μm^2 .< areas .< 500μm^2)
        else
            areas = lengths
            correct_size = trues(length(lengths))
            # remove background
            correct_size[1] = false
        end

        ids = slice_ids[correct_size]
        # get the centroids for all ids and then select only the ids that appear
        # in the current slice and *also* have the correct size
        centroids = Images.component_centroids(component_slice)[slice_ids[correct_size] .+ 1]

        n = length(centroids)
        ys = map(f->f[1], centroids) # y corresponds to rows
        xs = map(f->f[2], centroids) # x corresponds to columns
        frames = fill(idx-1, n)

        # select the current time slice and enforce storage order to match
        # components
        slice = view(img, Axis{:y}(:), Axis{:x}(:), Axis{:channel}(:), Axis{:time}(timepoint))
        localities = get_localities(component_slice, ids, dist=dist)
        # dictionary of ids to areas
        data = OrderedDict(:x=>xs,
                           :y=>ys,
                           :frame=>frames,
                           :id=>ids,
                           :area=>areas[correct_size])

        cax = AxisArrays.axes(img, Axis{:channel})
        for c in cax
            channelslice = view(slice, Axis{:channel}(c))
            tfs = Float64[]
            medbkgs = Float64[]

            for (id, indices) in localities
                # total fluorescence is the sum of all signal in the actual
                # footprint of the object. We have to do a copy operation here
                # due to https://github.com/JuliaArrays/AxisArrays.jl/issues/179
                push!(tfs, sum(channelslice[component_slice .== id]))

                # median background is the median of background signal in the
                # locality of object
                push!(medbkgs, median(channelslice[indices]))
            end
            data[Symbol("tf_", c)] = tfs
            data[Symbol("medbkg_", c)] = medbkgs
        end
        push!(particles, DataFrames.DataFrame(data))
    end
    vcat(particles...)
end

function build_tp_df(img::AxisArray{T1, 3},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Bool}
    build_tp_df(AxisArray(reshape(img, size(img)..., 1),
                          AxisArrays.axes(img)..., Axis{:channel}([:slice])),
                thresholds;
                dist=dist
               )
end

function build_tp_df(img::AxisArray{T1, 4},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Bool}

    # we have to pass the underlying array due to
    # https://github.com/JuliaImages/ImageMorphology.jl/issues/21
    components = Images.label_components(thresholds.data, [axisdim(thresholds, Axis{:y}), axisdim(thresholds, Axis{:x})])
    build_tp_df(img, AxisArray(components, AxisArrays.axes(thresholds)); dist=dist)
end

function build_tp_df(img::AxisArray{T1, 3},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Integer}

    build_tp_df(AxisArray(reshape(img, size(img)..., 1),
                          AxisArrays.axes(img)..., Axis{:channel}([:slice])),
                thresholds;
                dist=dist
               )
end

"""
Given an `img` with at least `y`, `x`, and `t` axes and a 3 dimensional boolean
array, `thresholds`, in yxt order.
"""
function build_tp_df(img::AxisArray{T1, 4},
                     thresholds::BitArray{3}; dist=(2, 10)) where {T1}
    _axes = Tuple(AxisArrays.axes(img, Axis{ax}) for ax in (:y, :x, :time))
    build_tp_df(img, AxisArray(Bool.(thresholds), _axes...), dist=dist)
end

"""
    get_locality(foreground, cell_mask; dist=(mindist, maxdist))

Given a boolean matrix of the pixels belonging to a single cell,
`cell_mask`, and a boolean matrix of foreground pixels, `foreground`, this
function  identifies the local background ring around the cell that is at
minimum `mindist` away from every object (in pixels) and a maximum of `maxdist`
away from the target cell.
"""
function get_locality(foreground::AbstractArray{Bool, 2}, cellmask::AbstractArray{Bool, 2}; dist=(2, 10))
    all_localities = dist[1] .< distance_transform(feature_transform(foreground)) .< dist[2]
    locality = dist[1] .< distance_transform(feature_transform(cellmask)) .< dist[2]
    # only return true where this cells locality overlaps with other localities, this way we
    # insure that only foreground areas are counted in the locality
    locality .&= all_localities
end

"""
    get_localities(labels, ids; dist) -> Dict

Gets the local background areas given the output of `label_components` and a
list of objects whose areas should be found defined by `ids`. The local
background is defined as a ring around the object that is at
minimum `mindist` away from every object (in pixels) and a maximum of `maxdist`
away from the target object. The result is a dictionary mapping the object id to
all indices in its local background.

!!! warning
    It's really important that `labels` has all foreground objects labeled
    because otherwise they might inadvertently be included in an object's
    locality.
"""
function get_localities(labels::AbstractArray{Int, 2}, ids::Vector{Int}; dist=(2, 10))
    boxes = component_boxes(labels)
    # boxes always contains the background label as well, while there's no guarantee
    # that ids does as well
    (!(0 in ids)) && deleteat!(boxes, 1)

    # if we detect gaps in the numbering it's probably because some labels were
    # removed, which can make our foreground detection incorrect
    computed_ids = sort(unique(labels))[2:end]
    if computed_ids != collect(minimum(computed_ids):maximum(computed_ids))
        @warn "Do not remove values from `labels`, only from `ids`. Calc might be wrong!"
    end

    nx, ny = size(labels)
    localities = OrderedDict{Int, Vector{CartesianIndex{2}}}()

    allobjects = labels .> 0

    δ = dist[1] + dist[2]

    for id in ids
        # unpack bounding box
        (minx, miny), (maxx, maxy) = boxes[id]

        # extend window around bounding box by the max distance
        xrange = max(minx-δ, 1):min(maxx+δ, nx)
        yrange = max(miny-δ, 1):min(maxy+δ, ny)

        local_allobjects = view(allobjects, xrange, yrange)
        local_cellmask = view(labels, xrange, yrange) .== id
        @assert sum(local_cellmask) > 0 "No object of $id found"

        locality = get_locality(local_allobjects, local_cellmask; dist=dist)

        # use offset arrays to return CartesianIndices in the original image coordinates
        localities[id] = findall(OffsetArray(locality, xrange, yrange))
    end

    localities
end