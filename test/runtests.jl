using Test
using Distributions
using Colors
using FixedPointNumbers
using SegmentationTools
using ImageDraw
using Images
using AxisArrays
using Unitful: μm

@testset "Find peak centers" begin
    # even width peak
    xs = Float64.(1:100)
    ys = fill(0.0, 100)
    ys[20:29] .+= 2
    @test SegmentationTools.find_peak_center(xs, ys) == 24.5

    # odd width peak
    xs = Float64.(1:100)
    ys = fill(0.0, 100)
    ys[20:30] .+= 2
    @test SegmentationTools.find_peak_center(xs, ys) == 25.0

    # double peaks
    xs = Float64.(1:100)
    ys = fill(0.0, 100)
    ys[20:30] .+= 2
    ys[50:60] .+= 2
    @test_throws ErrorException("Dependent variable is not unimodal") SegmentationTools.find_peak_center(xs, ys) == 25.0

    # shoulder
    xs = Float64.(1:100)
    ys = fill(0.0, 100)
    ys[20:30] .+= 1
    ys[26:30] .+= 1
    # high cut
    @test SegmentationTools.find_peak_center(xs, ys) == 28.0
    # low cut
    @test SegmentationTools.find_peak_center(xs, ys, h=0.25) == 25.0
end

@testset "colorize" begin
    # Fixed Point Number
    data = rand(Gray{N0f16}, 1024, 1024, 10)
    # make sure the backing type is still N0f16 changed by the mapping
    @test eltype(SegmentationTools.colorize(data)) <: RGB{N0f16}

    # same for floating point numbers
    data = rand(Gray{Float64}, 1024, 1024, 10)
    @test eltype(SegmentationTools.colorize(data)) <: RGB{Float64}
end

@testset "trackpy table construction" begin
    img = zeros(Gray, 100, 100, 5)

    # create three cells
    centers = [Point(20, 20), Point(25, 17), Point(70, 70), Point(70, 70), Point(16, 80)]
    radii = [5, 5, 15, 5, 7]

    bkg = 0.1

    for i in 1:size(img, 3)
        # draw the "cells" and make them grow dimmer over time
        colors = [Gray(1.0/i) for _ in centers]
        # the 3rd circle is concentric with the 3rd and is a high background
        colors[3] = Gray(bkg)
        ImageDraw.draw!(view(img, :, :, i), [ImageDraw.CirclePointRadius(c, r) for (c,r) in zip(centers, radii)], colors)
    end

    wrapped = AxisArray(img, Axis{:y}(1μm:1μm:100μm), Axis{:x}(1μm:1μm:100μm), Axis{:time}(1:5));

    mask = AxisArray(img .> bkg, AxisArrays.axes(wrapped))
    particles = SegmentationTools.build_tp_df(wrapped, mask)

    # three timepoints and three separate particles
    @test size(particles, 1) == 3*5

    # we now assign unique ids per-frame, so we need to convert these back to
    # consistent ids across time
    ids = mod.(particles[!, :id].-1, 3) .+ 1
    # verify that the total fluorescence values in the dataframe grow dimmer as
    # we expect
    @test all(particles[ids .== 1, :tf_slice] ./ 145.0 .≈ [1.0/i for i in 1:5])

    # for the third particle verify that the total fluorescence values have the
    # correct amount of background subtracted off, i.e for a background value of
    # `bkg` the final value should be tf - area * bkg
    @test all(particles[ids .== 3, :medbkg_slice] .≈ [bkg for i in 1:5])

    # test centroids
    @test particles[1, :x] == centers[5].x
    @test particles[1, :y] == centers[5].y
end

@testset "local background calc" begin
    include("localities.jl")
end

@testset "Flatfield correction" begin
    # construct a 2D gaussian
    D = reshape(pdf.(Normal(50, 100), range(0, length=100)), :, 1) .* 100
    # 0-centered gaussian noise
    noise = Gray.(rand(Normal(0, 0.001), 100, 100))

    illum = Gray.(D * D')
    # make the signal = 2*illumination intensity
    foreground = AxisArray(illum .+ illum .+ noise)
    # background
    background = illum .+ noise

    # if the flatfield correction is working properly than the image after
    # correction should equal 2.0 everywhere since it's 2 * illum
    @test all(SegmentationTools.flatfield_correct(foreground, background, noise) .≈ 2.0)
end