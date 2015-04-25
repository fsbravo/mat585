%%% calculate affinity between each pair of images using a gaussian kernel

function [A] = gaussian_kernel_weights(imgs)

    n = size(imgs,3);
    A = zeros(n);
    for i=1:n
        img_i = imgs(:,:,i);
        for j=i:n
            img_j = imgs(:,:,j);
            A(i,j) = exp(-norm(img_i(:)-img_j(:))^2);
        end
    end
    A = A + triu(A,1)';
    
end