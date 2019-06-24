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
