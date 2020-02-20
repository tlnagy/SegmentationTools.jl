img = zeros(Gray, 50, 50);

# create a cell with a background half a bright around it
centers = [Point(20, 20), Point(20, 20)]
radii = [15, 5]
# the local background will be twice as dim as the "cell"
colors = [Gray(0.5), Gray(1.0)]

ImageDraw.draw!(img, [ImageDraw.CirclePointRadius(c, r) for (c,r) in zip(centers, radii)], colors)
labels = Images.label_components(img .== 1.0)
localities = SegmentationTools.get_localities(labels, [1], dist=(0, 10))

@test median(Float64.(img[localities[1]])) == 0.5

img = zeros(Gray, 50, 50);
img[10, 10] = Gray(1.0);
img[20, 10] = Gray(1.0);
img[30, 10] = Gray(1.0);

labels = Images.label_components(img .== 1.0)
labels[labels .== 2] .= 0
ids = [1, 3]

msg = "Do not remove values from `labels`, only from `ids`. Calc might be wrong!"
@test_logs (:warn, msg) get_localities(labels, ids; dist=(5, 15))

labels = Images.label_components(img .== 1.0)

ids = [1]
localities = get_localities(labels, ids, dist=(5,50))

@test all(keys(localities) .== ids)

locality = Set(localities[1])

# all of the points should be missing from the locality
@test length(intersect(findall(img .> 0), locality)) == 0