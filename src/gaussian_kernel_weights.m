%%% calculate affinity between each pair of images using a gaussian kernel

function [W] = gaussian_kernel_weights(imgs)
    n = size(imgs,3);
    W = zeros(n);
    for i=1:n
        img_i = imgs(:,:,i);
        for j=i:n
            img_j = imgs(:,:,j);
            distances(i,j) = norm(img_i(:)-img_j(:))^2;  %this is on the order of 1e6
        end
    end
    distances = distances + triu(distances,1)';
    W = exp(-distances./median(distances(:)));
end