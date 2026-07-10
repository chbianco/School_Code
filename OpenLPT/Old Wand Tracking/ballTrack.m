im = imread('wand_images/cam1/frame_00200.tif');
[centers, radii, metric] = imfindcircles(im, [20, 90], 'ObjectPolarity', 'dark');

metricStrong2 = metric(1:2)
centersStrong2 = centers(1:2,:)
radiiStrong2 = radii(1:2)

imshow(im)
viscircles(centersStrong2, radiiStrong2,'EdgeColor','b');
