function [image_matrix] = image_reader(set_name,range)
% image_reader.m
% parse tif images from one of the three datasets.
% INPUT:
%   set_name: either drosophilia_fixed, drosophilia_live, or zebrafish
%   The relative path to these directories should be ../images/set_name
%   range: integer = the last filename number. 

% OUTPUT:
%   image_matrix: 100x100xrange 3D array containing ordered, registered
%   images

switch set_name
    case 'drosophilia_fixed'
        path = '../images/drosophilia_fixed/';
    case 'drosophilia_live'
        path = '../images/drosophilia_live/';
    case 'zebrafish'
        path = '../images/zebrafish/';
    otherwise
        error('No such data set! Use drosophilia_fixed, drosophilia_live, or zebrafish');
end

for i=1:range
    if i < 10
       filename = strcat('image0',int2str(i),'.tif'); 
    else
       filename = strcat('image',int2str(i),'.tif');   
    end
    filepath = strcat(path,filename); 
    %should probably preallocate, but we need to know the size of the
    %images beforehand
    image_matrix(:,:,i) = imread(filepath,'tiff');
end
% grayscale images should be MxN matrix
%  For TIFF files containing color images that use the CMYK color space,
%  A is an M-by-N-by-4 array. 
end