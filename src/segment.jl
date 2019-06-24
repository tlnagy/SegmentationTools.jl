using Images: otsu_threshold, imadjustintensity

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
    build_tp_df(img, segments, signal_channel)

Build a `DataFrames.DataFrame` that is compatible with trackpys `link_df` function. Needs
to be converted to a `Pandas.DataFrame` before passing
"""
function build_tp_df(img::AbstractArray{T, 4},
                     segments::Array{S, 1},
                     signal_channel) where {T, S <: SegmentedImage}

    particles = DataFrames.DataFrame[]
    for i in 1:length(segments)
        # drop first element, which is the background
        centroids =  component_centroids(segments[i].image_indexmap)[2:end]
        n = length(centroids)
        ys = map(f->f[1], centroids) # y corresponds to rows
        xs = map(f->f[2], centroids) # x corresponds to columns
        frames = fill(i-1, n)

        signal = img[Axis{:channel}(signal_channel), Axis{:time}(i)]
        tf = fill(0.0, n)
        for obj in 1:n
            tf[obj] = sum(signal[segments[i].image_indexmap .== obj])
        end
        ids = segments[i].segment_labels
        # dictionary of ids to areas
        areas = segments[i].segment_pixel_count
        data = Dict(:x=>xs,
                    :y=>ys,
                    :frame=>frames,
                    :id=>ids,
                    :area=>map(id->areas[id], ids),
                    :tf=>tf)
        push!(particles, DataFrames.DataFrame(data))
    end
    vcat(particles...)
end
