"""
    flatfield(img)

Compute the flatfield image from `img` by calculating the median value
along the 3rd dimension which normally corresponds to different
positions. This way the contribution of the optical path can be
determined. Requires an `img` with labeled dimensions `x` and `y`, order
does not matter.
"""
function flatfield(img::AxisArray{T, 3}) where {T}
    ysize, xsize = size(img, Axis{:y}), size(img, Axis{:x})
    med = Array{T}(undef, ysize, xsize)

    R = CartesianIndices((xsize, ysize))
    for i in R
        med[i] = median(img[Axis{:y}(i[1]), Axis{:x}(i[2])])
    end
    med
end

"""
	flatfield_correct(img, axis, darkfield, flatfield)

Flatfield correct `img`. This function subtracts the darkfield, clamps the output
to positive values, divides by the flatfield, and then rescales and rebuilds the
original image perserving the properties and axes.

For example, to correct along the DAPI channel:

```
flatfield_correct(img[Axis{:channel}(:EPI_DAPI)], darkfield, flatfield_dapi)
```
"""
function flatfield_correct(img::AxisArray{T1, N},
                           darkfield::AbstractArray{T2, 2},
                           flatfield::AbstractArray{T3, 2}) where {T1, T2, T3, N}
    out = T1.(imadjustintensity(clamp01.(Float32.(arraydata(img) .- darkfield))./flatfield))
    AxisArray(out, img.axes)
end

function flatfield_correct(img::ImageMeta{T, N},
                           darkfield,
                           flatfield) where {T, N}
    out = flatfield_correct(img.data, darkfield, flatfield)
    copyproperties(img, out)
end

light_source_contrib(img::AxisArray, seeds::BitArray{3}) = light_source_contrib(img, AxisArray(seeds, img.axes))

"""
    light_source_contrib(img, seeds, dist, h)

Computes the contribution of the light source fluctuations for each time point by first
identifying background pixels by finding pixels that are more than `dist` away from true
values in `seed`. A kernel density estimate is then fit to the background pixels to get
a continuous distribution and then [`find_peak_center`](@ref) for more details how the "center" is defined. 

### Rationale

Certain light sources, especially arc lamps, exhibit substantial variance in 
total brightness in time. Additionally, the background signal is dependent on the total
light delivered and we can use that to identify the arc lamp wander. 


### Example:

```julia
a = img[Axis{:position}(2), Axis{:channel}(:EPI_mNG)].data
b = img[Axis{:position}(2), Axis{:channel}(:EPI_BFP)].data .> 0.01
SegmentationTools.light_source_contrib(a, b)
```
"""
function light_source_contrib(img::AxisArray{T1, 3},
                              seeds::AxisArray{T2, 3};
                              dist::Int = 30) where {T1, T2 <: Bool}
    n = length(timeaxis(img))
    (n != length(timeaxis(seeds))) && error("The time dimensions of both arrays need to the same")
    bkg_means = Array{Float64}(undef, n)
    for i in 1:n
        seed_slice = seeds[Axis{:time}(i)]
        bkg_mask = distance_transform(feature_transform(seed_slice)) .> dist
        # fit a kernel density estimate to the background pixels
        fitted_kde = kde(vec(Float64.(img[Axis{:time}(i)][bkg_mask])))
        # find the center of the KDE of the background peak
        bkg_means[i] = find_peak_center(collect(fitted_kde.x), fitted_kde.density)
    end
    bkg_means
end

"""
    find_peak_center(x, y, h)

Given two equal length arrays, `x`, and `y`, computes the half width 
at `h` relative to the maximum of the peak in `y`, i.e. for `h=0.5`, this 
function returns the value of `x` that corresponds to the half-width 
half-max. This is a robust method for identifying the center of a peak.
"""
function find_peak_center(x::Vector{Float64}, y::Vector{Float64}; h::Float64=0.75)
    (length(x) != length(y)) && throw(DimensionMismatch("x and y have to have equal lengths"))
    # value at which we cut the peak
    cutoff = maximum(y)*h
    
    # find the indices where we "enter" and "exit" the peak
    transitions = findall(i-> i != 0, diff(y .> cutoff))
    (length(transitions) != 2) && error("Dependent variable is not unimodal")

    # translate from indices to the original units of `slice`
    low = x[transitions[1] + 1]
    high = x[transitions[2]]
    low + (high - low) / 2
end