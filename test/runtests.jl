using Test
using Distributions
using SegmentationTools

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
