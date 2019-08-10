var documenterSearchIndex = {"docs":
[{"location":"#Segment-cells-with-Julia-1","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"","category":"section"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"Welcome to SegmentationTools.jl! This is a random collection of functions that I've implemented for segmenting single cells from microscopy images.","category":"page"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"CurrentModule=SegmentationTools","category":"page"},{"location":"#Flow-1","page":"Segment cells with Julia","title":"Flow","text":"","category":"section"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"The general procedure follows the following format:","category":"page"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"Segmentation -> Linking -> Analysis","category":"page"},{"location":"#Segmentation-1","page":"Segment cells with Julia","title":"Segmentation","text":"","category":"section"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"SegmentationTools.segment_cells","category":"page"},{"location":"#SegmentationTools.segment_cells","page":"Segment cells with Julia","title":"SegmentationTools.segment_cells","text":"segment_cells(img, seed_channel, segment_channel, maskfunc)\n\n\n\n\n\n","category":"function"},{"location":"#Linking-1","page":"Segment cells with Julia","title":"Linking","text":"","category":"section"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"SegmentationTools.build_tp_df","category":"page"},{"location":"#SegmentationTools.build_tp_df","page":"Segment cells with Julia","title":"SegmentationTools.build_tp_df","text":"build_tp_df(img, segments, signal_channel)\n\nBuild a DataFrames.DataFrame that is compatible with trackpys link_df function. Needs to be converted to a Pandas.DataFrame before passing\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"I generally use it like so:","category":"page"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"particles = SegmentationTools.build_tp_df(img1, segments, :EPI_mNG);\nusing PyCall\ntp = pyimport(\"trackpy\")\n\nimport Pandas\n# convert particles to a Pandas DataFrame and call trackpy\nt = tp.link_df(Pandas.DataFrame(particles), 5);\nlinked = DataFrame(Pandas.DataFrame(t)); # convert back to Julian dataframe","category":"page"},{"location":"#Analysis-and-Plotting-1","page":"Segment cells with Julia","title":"Analysis and Plotting","text":"","category":"section"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"SegmentationTools.create_cell_grid","category":"page"},{"location":"#SegmentationTools.create_cell_grid","page":"Segment cells with Julia","title":"SegmentationTools.create_cell_grid","text":"create_cell_grid(img, tracks; win) -> grid\n\nGiven a xyt image and set of tracks, creates a grid of all tracked cells with the cells centered in their grid location. Helps to quickly diagnose weird cell segmentation behavior.\n\n\n\n\n\n","category":"function"},{"location":"#Miscellaneous-1","page":"Segment cells with Julia","title":"Miscellaneous","text":"","category":"section"},{"location":"#","page":"Segment cells with Julia","title":"Segment cells with Julia","text":"SegmentationTools.find_peak_center\nSegmentationTools.light_source_contrib\nSegmentationTools.flatfield_correct","category":"page"},{"location":"#SegmentationTools.find_peak_center","page":"Segment cells with Julia","title":"SegmentationTools.find_peak_center","text":"find_peak_center(x, y, h)\n\nGiven two equal length arrays, x, and y, computes the half width  at h relative to the maximum of the peak in y, i.e. for h=0.5, this  function returns the value of x that corresponds to the half-width  half-max. This is a robust method for identifying the center of a peak.\n\n\n\n\n\n","category":"function"},{"location":"#SegmentationTools.light_source_contrib","page":"Segment cells with Julia","title":"SegmentationTools.light_source_contrib","text":"light_source_contrib(img, seeds, dist, h)\n\nComputes the contribution of the light source fluctuations for each time point by first identifying background pixels by finding pixels that are more than dist away from true values in seeds. A kernel density estimate is then fit to the background pixels to get a continuous distribution and then find_peak_center is used to identify the \"average\" response of a background pixel.\n\nnote: Note\nThe distribution of background pixels often ends up being quite non-normal, likely due to misidentification. Empirically, I've found the finding the center of the kernel  density estimate to be more robust than simpler statistics like mean, median, or even the maximum value of the KDE.\n\nRationale\n\nCertain light sources, especially arc lamps, exhibit substantial variance in  total brightness in time. Additionally, the background signal is dependent on the total light delivered and we can use that to identify the arc lamp wander and subtract it from the whole field of view at each timestep. This gives us more stable total fluorescence values over time.\n\nExample:\n\na = img[Axis{:position}(2), Axis{:channel}(:EPI_mNG)].data\nb = img[Axis{:position}(2), Axis{:channel}(:EPI_BFP)].data .> 0.01\nSegmentationTools.light_source_contrib(a, b)\n\n\n\n\n\n","category":"function"},{"location":"#SegmentationTools.flatfield_correct","page":"Segment cells with Julia","title":"SegmentationTools.flatfield_correct","text":"flatfield_correct(img, axis, darkfield, flatfield)\n\nFlatfield correct img. This function subtracts the darkfield, clamps the output to positive values, divides by the flatfield, and then rescales and rebuilds the original image perserving the properties and axes.\n\nFor example, to correct along the DAPI channel:\n\nflatfield_correct(img[Axis{:channel}(:EPI_DAPI)], darkfield, flatfield_dapi)\n\n\n\n\n\n","category":"function"}]
}
