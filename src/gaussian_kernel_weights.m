%%% calculate affinity between each pair of images using a gaussian kernel

function [A] = gaussian_kernel_weights(imgs)
    n = size(imgs,3);
    A = zeros(n);
    for i=1:n
        img_i = imgs(:,:,i);
        for j=i:n
            img_j = imgs(:,:,j);
            difference = norm(img_i(:)-img_j(:))^2;  %this is on the order of 1e6
            %normalize by average of squared frobenius norms
            average_norm = (norm(img_i(:))^2 + norm(img_j(:))^2)/2;
            A(i,j) = exp(-difference/average_norm);
        end
    end   
    A = A + triu(A,1)';
end