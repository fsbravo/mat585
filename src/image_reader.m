function [image_mat, nchannels] = image_reader(set_name,range,varargin)
% image_reader.m
% parse tif images from one of the three datasets.
% INPUT:
%   set_name: either drosophilia_fixed, drosophilia_live, or zebrafish
%   The relative path to these directories should be ../images/set_name
%   range: integer = the last filename number. 

% OUTPUT:
%   image_matrix: 100x100xrange 3D array containing ordered, registered
%   images

% grayscale images should be MxN matrix
%  For TIFF files containing color images that use the CMYK color space,
%  A is an M-by-N-by-4 array. 
%  For TIFF files containing color images that use the RGB color space,
%  A is an M-by-N-by-4 array. 

if ~isempty(varargin)
    indicator_rgb = varargin{1};
else
    indicator_rgb = 0;
end
% indicator_rgb = 0; %for drosophila_fixed, 0 for grayscale, 1 for RGB
nchannels = 1; %grayscale

switch set_name
    case 'drosophila_fixed'
        path = '../images/drosophila_fixed/';
        %These images have 3 RGB channels
        if (indicator_rgb)
            nchannels = 3;
        end
    case 'drosophila_live'
        path = '../images/drosophila_live/';
    case 'zebrafish'
        path = '../images/zebrafish/';
    otherwise
        error('No such data set! Use drosophila_fixed, drosophila_live, or zebrafish');
end
%We need to know the dimensions of the images a priori to preallocate
image_matrix = zeros(100,100,nchannels,range);
for i=1:range
    if i < 10
       filename = strcat('image0',int2str(i),'.tif'); 
    else
       filename = strcat('image',int2str(i),'.tif');   
    end
    filepath = strcat(path,filename); 
    if (strcmp(set_name,'drosophila_fixed') && indicator_rgb == 0)
        image_matrix(:,:,1,i) = rgb2gray(imread(filepath,'tiff'));
    else
        image_matrix(:,:,:,i) = imread(filepath,'tiff');
    end
end

%%% preprocess images

switch set_name
    case 'drosophila_fixed' 
        if (indicator_rgb ~= 1)
            error('cant preprocess grayscale drosophila images')
        end
        channel_weight = [0.5 1 1];
        channel_blur = [0.05 0.05 0.05];
        channel_normalize = [1 0 0];
        channel_mean_center = [1 0 0];
        resize_image = true;
    case 'drosophila_live'
        channel_weight = 1;
        channel_blur = 0.05;
        channel_normalize = 1;
        channel_mean_center = 0;
        resize_image = false;
    case 'zebrafish'
        channel_weight = 1;
        channel_blur = 0;
        channel_normalize = 0;
        channel_mean_center = 1;
        resize_image = false;
    otherwise
        error('No such data set! Use drosophila_fixed, drosophila_live, or zebrafish');
end
image_mat = apply_image_functions(image_matrix, 100, 2, channel_weight, ...
    channel_blur, channel_normalize, channel_mean_center, resize_image);
end