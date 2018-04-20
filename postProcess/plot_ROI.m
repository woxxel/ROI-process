


ROIs = h5read('/media/mizuta/Analyze_AS1/245/Session01/imageStack/ImagingData_main.h5','/MERGED/A');

imshow(sum(ROIs,3));