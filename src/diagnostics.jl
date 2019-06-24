"""
r   create_cell_grid(img, tracks; win) -> grid

Given a xyt image and set of tracks, creates a grid of all tracked cells with the cells centered
in their grid location. Helps to quickly diagnose weird cell segmentation behavior.


"""
function create_cell_grid(img::AbstractArray{T, 3}, tracks::AbstractArray{S, 1}; win=30) where {T, S <: SubDataFrame}
    n_tracks = length(tracks)
    ny, nx, nt = size(img)
    cols = ceil(Int, sqrt(n_tracks))
    fontface = newface(normpath(@__DIR__, "..", "assets", "ubuntu-font-family-0.83", "Ubuntu-R.ttf"))
    w = cols*(win*2+2)
    h = ceil(Int, n_tracks / cols)*(win*2+2)
    # the default element should be different for binary masks vs greyscales
    # this picks something intelligent for both
    default_el = T <: Colorant ? maximum(img[:, :, 1]) : oneunit(UInt8)
    grid = fill(default_el, h, w, nt)

    for (pidx, subdf) in enumerate(tracks)
        # calculate the x,y center of this grid location
        xcent, ycent = (pidx-1)÷cols*(win*2+2)+win+1, mod((pidx-1)*(win*2+2)+win+1, w)
        particle = subdf[1, :particle]

        # for each frame in the movie
        for (tidx, xpos, ypos, id) in zip(1:nt, subdf[:x], subdf[:y], subdf[:id])

            # find approximate index of the center of the cell
            xidx, yidx = round(Int, ypos), round(Int, xpos)
            # handle edges intelligently
            xrange, yrange = max(xidx-win, 1):min(xidx+win, nx), max(yidx-win, 1):min(yidx+win, ny)

            # range - idx centers the view on the cell to be +/- the window and then we offset it by
            # our current location in the grid
            xloc, yloc = xrange .- xidx .+ xcent,yrange .- yidx .+ ycent

            if T <: Colorant
                cell = img[xrange, yrange, tidx]
            else # if a color wasn't provided then we need to select only cells that are labeled
                cell = UInt8.(img[xrange, yrange, tidx] .== id)
            end

            renderstring!(cell, "$particle", fontface, (10,10), 2, 2, halign=:hleft, valign=:vtop, fcolor=default_el)

            grid[xloc, yloc, tidx] .= cell
        end
    end
    grid
end
