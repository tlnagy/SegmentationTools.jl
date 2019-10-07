using Images: otsu_threshold, imadjustintensity
using Unitful: μm
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
    build_tp_df(img, thresholds; dist)

Build a `DataFrames.DataFrame` that is compatible with trackpys `link_df`
function. Needs to be converted to a `Pandas.DataFrame` before passing to
trackpy. `dist` is a 2-tuple of integers indicating the minimum and maximum
distance away in pixels from each cell to include in its local background
calculation.
"""
function build_tp_df(img::AxisArray{T1, 4},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Bool}

    xstep = step(AxisArrays.axes(img, Axis{:x}).val)
    ystep = step(AxisArrays.axes(img, Axis{:y}).val)
    (xstep != ystep) && @warn "Different scaling for x and y axes is not supported"
    pixelarea = xstep * ystep
    particles = DataFrames.DataFrame[]
    # we have to pass the underlying array due to
    # https://github.com/JuliaImages/ImageMorphology.jl/issues/21
    components = Images.label_components(thresholds[Axis{:y}(:),
                                                    Axis{:x}(:),
                                                    Axis{:time}(:)].data, [1,2])
    @showprogress 1 "Computing..." for (idx, timepoint) in enumerate(timeaxis(img))
        component_slice = view(components, :, :, idx)

        # convert boolean area to microns with the assumption that the pixel
        # space is the same in both x and y
        lengths = Images.component_lengths(component_slice)
        areas = round.(μm^2, lengths .* pixelarea, sigdigits=4)
        # filter out too large and too small particles
        correct_size = (10μm^2 .< areas .< 500μm^2)
        ids = unique(component_slice)[correct_size[lengths .> 0]]
        centroids = Images.component_centroids(component_slice)[correct_size]

        n = length(centroids)
        ys = map(f->f[1], centroids) # y corresponds to rows
        xs = map(f->f[2], centroids) # x corresponds to columns
        frames = fill(idx-1, n)

        # select the current time slice and enforce storage order to match
        # components
        slice = view(img, Axis{:y}(:), Axis{:x}(:), Axis{:channel}(:), Axis{:time}(timepoint))
        tf, bkg = _compute_total_fluorescences(slice, ids, component_slice, centroids, dist=dist)

        # dictionary of ids to areas
        data = Dict(:x=>xs,
                    :y=>ys,
                    :frame=>frames,
                    :id=>ids,
                    :area=>areas[correct_size])

        cax = AxisArrays.axes(img, Axis{:channel})
        for c in 1:size(tf, 2)
            data[Symbol("tf_", cax[c])] = tf[:, c]
            data[Symbol("bkg_", cax[c])] = bkg[:, c]
            data[Symbol("normed_tf_", cax[c])] = tf[:, c] .- bkg[:, c]
        end
        push!(particles, DataFrames.DataFrame(data))
    end
    vcat(particles...)
end

function build_tp_df(img::AxisArray{T1, 3},
                     thresholds::AxisArray{T2, 3}; dist=(2, 10)) where {T1, T2 <: Bool}
    build_tp_df(AxisArray(reshape(img, size(img)..., 1),
                          AxisArrays.axes(img)..., Axis{:channel}([:slice])),
                thresholds,
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
    _get_locality_mask(cell_mask, foreground; dist=(mindist, maxdist))

    Given a boolean matrix of the pixels belonging to a single cell,
`cell_mask`, and a boolean matrix of foreground pixels, `foreground`, this
function  identifies the local background ring around the cell that is at
minimum `mindist` away from every cell and a maximum of `maxdist` away from the
target cell.
"""
function _get_locality_mask(cell_mask::BitArray{2}, foreground::BitArray{2}; dist=(1,6))
    all_localities = dist[1] .< distance_transform(feature_transform(foreground)) .< dist[2]
    # create a 5px wide mask that is at least 1 px away from the cell
    locality = dist[1] .< distance_transform(feature_transform(cell_mask)) .< dist[2]
    # only return true where this cells locality overlaps with all local areas
    locality .&= all_localities
end


"""
    _compute_equivalent_background(slice, labels, id; dist)

    Given a cell `id` and a matrix of labeled cells `labels`, compute the total
fluorescence expected for an object the size of the cell using the local background
fluorescence surrounding the cell. The local background is computed from `dist[1]`
to `dist[2]` away from the cell.
"""
function _compute_equivalent_background(slice::AbstractArray{T, 2},
                                       labels::AbstractArray{Int, 2},
                                       id::Int;
                                       dist=(2, 10)) where {T}
    foreground = labels .!= 0.0
    component_mask = labels .== id
    locality_mask = _get_locality_mask(component_mask, foreground; dist=dist)
    local_bkg = locality_mask .* slice
    # fluorescence of an equivalent background area
    median(local_bkg[locality_mask]) * count(component_mask)
end

function _compute_total_fluorescences(slice::AxisArray{T1, 3},
                                      ids::Vector{Int},
                                      labels::AbstractArray{Int, 2},
                                      centroids::AbstractArray{Tuple{Float64, Float64}};
                                      dist=(2, 10)) where {T1}
    n = length(centroids)
    nₛ = size(slice, Axis{:channel})
    tf = fill(0.0, n, nₛ)
    bkg = fill(0.0, n, nₛ)
    for c in 1:nₛ
        signal = view(slice, Axis{:channel}(c))
        _tf, _bkg = _compute_total_fluorescences(signal, ids, labels, centroids, dist=dist)
        tf[:, c] .= _tf
        bkg[:, c] .= _bkg
    end
    tf, bkg
end

function _compute_total_fluorescences(slice::AxisArray{T1, 2},
                                      ids::Vector{Int},
                                      labels::AbstractArray{Int, 2},
                                      centroids::AbstractArray{Tuple{Float64, Float64}};
                                      dist=(2, 10)) where {T1}
    n = length(centroids)
    nx = size(slice, Axis{:x})
    ny = size(slice, Axis{:y})
    win = 50
    tf = fill(0.0, n)
    bkg = fill(0.0, n)

    for (idx, id) in enumerate(ids)
        # compute and subtract the local bkg equivalent from the total
        # fluorescence
        centroid = centroids[idx]
        xidx, yidx = round(Int, centroid[2]), round(Int, centroid[1])
        xrange = max(xidx-win, 1):min(xidx+win, nx)
        yrange = max(yidx-win, 1):min(yidx+win, ny)
        local_signal = view(slice, Axis{:y}(yrange), Axis{:x}(xrange))
        local_labels = view(labels, yrange, xrange)
        bkg[idx] = _compute_equivalent_background(local_signal, local_labels, id, dist=dist)
        tf[idx] = sum(slice[labels .== id])
    end

    tf, bkg
end