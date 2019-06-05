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
