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
    radii = [5, 5, 11, 5, 7]

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
    # verify that the total fluorescence values in the dataframe grow dimmer as
    # we expect
    @test all(particles[particles[!, :id] .== 1, :tf_slice] ./ 145.0 .≈ [1.0/i for i in 1:5])

    # for the third particle verify that the total fluorescence values have the
    # correct amount of background subtracted off, i.e for a background value of
    # `bkg` the final value should be tf - area * bkg
    @test all(particles[particles[!, :id] .== 3, :tf_slice] .≈ [69/i - 69*bkg for i in 1:5])

    # test centroids
    @test particles[1, :x] == centers[5].x
    @test particles[1, :y] == centers[5].y

    # Check warnings
    wrapped = AxisArray(img, Axis{:y}(0.5μm:0.5μm:50μm), Axis{:x}(1μm:1μm:100μm), Axis{:time}(1:5));
    @test_logs (:warn, "Different scaling for x and y axes is not supported") SegmentationTools.build_tp_df(wrapped, mask)
end

@testset "local background calc" begin
    img = zeros(Gray, 50, 50)

    # create a cell with a background half a bright around it
    centers = [Point(20, 20), Point(20, 20)]
    radii = [15, 5]
    # the local background will be twice as dim as the "cell"
    colors = [Gray(0.5), Gray(1.0)]

    ImageDraw.draw!(img, [ImageDraw.CirclePointRadius(c, r) for (c,r) in zip(centers, radii)], colors)
    labels = Images.label_components(img .== 1.0)
    # the cell should be twice as bright as an equivalent sized background
    cell_tf = sum(img[img .== 1.0])
    @test cell_tf/SegmentationTools._compute_equivalent_background(img, labels , 1) == 2.0
end