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

@traitfn function flatfield_correct!(output,
                                     img::AA,
                                     flatfield::AbstractArray{T1, 2},
                                     darkfield::AbstractArray{T2, 2}) where {AA <: AxisArray, T1, T2; HasTimeAxis{AA}}
    for timepoint in timeaxis(img)
        tp = Axis{:time}(timepoint)
        _img = @view img[tp]
        _output = @view output[tp]
        flatfield_correct!(_output, _img, flatfield, darkfield)
    end
end

@traitfn function flatfield_correct!(output,
                                     img::AA,
                                     flatfield::AbstractArray{T1, 2},
                                     darkfield::AbstractArray{T2, 2}) where {AA <: AxisArray, T1, T2; !HasTimeAxis{AA}}
    R = CartesianIndices(Base.axes(img))
    for I in R
        output[I] = (float(img[I]) - darkfield[I]) / (flatfield[I] - darkfield[I])
    end
end

"""
    flatfield_correct(img, flatfield, darkfield)

Computes the flatfield for the given image `img` by dividing the image by the
flatfield after subtracting the darkfield from both[^1]

!!! note

    This function assumes that the flatfield image has *not* had the darkfield
    image subtracted yet.

[^1]: http://nic.ucsf.edu/resources/how-to-acquire-flat-field-correction-images/
"""
function flatfield_correct(img::AxisArray, flatfield::AbstractArray{T1, 2}, darkfield::AbstractArray{T2, 2}) where {T1, T2}
    output = similar(img, Gray{Float64})
    flatfield_correct!(output, img, flatfield, darkfield)
    output
end

flatfield_correct(img::ImageMeta, flatfield, darkfield) = copyproperties(img, flatfield_correct(img.data, flatfield, darkfield))


light_source_contrib(img::AxisArray, seeds::BitArray{3}; dist::Int=30, h::Float64=0.9) =
    light_source_contrib(img, AxisArray(seeds, img.axes); dist=dist, h=h)

"""
    light_source_contrib(img, seeds, dist, h)

Computes the contribution of the light source fluctuations for each time point by first
identifying background pixels by finding pixels that are more than `dist` away from true
values in `seeds`. A kernel density estimate is then fit to the background pixels to get
a continuous distribution and then [`find_peak_center`](@ref) is used to identify the
"average" response of a background pixel.

!!! note

    The distribution of background pixels often ends up being quite non-normal, likely due
    to misidentification. Empirically, I've found the finding the center of the kernel
    density estimate to be more robust than simpler statistics like mean, median, or even
    the maximum value of the KDE.

### Rationale

Certain light sources, especially arc lamps, exhibit substantial variance in
total brightness in time. Additionally, the background signal is dependent on the total
light delivered and we can use that to identify the arc lamp wander and subtract it
from the whole field of view at each timestep. This gives us more stable total
fluorescence values over time.


### Example:

```julia
a = img[Axis{:position}(2), Axis{:channel}(:EPI_mNG)].data
b = img[Axis{:position}(2), Axis{:channel}(:EPI_BFP)].data .> 0.01
SegmentationTools.light_source_contrib(a, b)
```
"""
function light_source_contrib(img::AxisArray{T1, 3},
                              seeds::AxisArray{T2, 3};
                              dist::Int = 30,
                              h::Float64 = 0.9) where {T1, T2 <: Bool}
    n = length(timeaxis(img))
    (n != length(timeaxis(seeds))) && error("The time dimensions of both arrays need to the same")
    bkg_means = Array{Float64}(undef, n)
    for i in 1:n
        seed_slice = seeds[Axis{:time}(i)]
        bkg_mask = distance_transform(feature_transform(seed_slice)) .> dist
        # fit a kernel density estimate to the background pixels
        fitted_kde = kde(vec(Float64.(img[Axis{:time}(i)][bkg_mask])))
        # find the center of the KDE of the background peak
        bkg_means[i] = find_peak_center(collect(fitted_kde.x), fitted_kde.density; h=h)
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


"""
    colorize(img, scheme)

Given a grayscale image `img`, apply a colorscheme `scheme` to the image lazily.
"""
function colorize(img::AbstractArray{Gray{T}, N}; scheme::Symbol=:magma) where {T, N}
    minval, maxval = gray.(extrema(img))
    cscheme = getfield(ColorSchemes, scheme)
    mappedarray(x->RGB{T}(get(cscheme, gray(x), (minval, maxval))), img)
end