function [diff_map] = full_diffusion_map(W,t)
    n = size(W,1);
    D_inv = diag(1./sum(W,2));
    A = D_inv*W;
    [V,E] = eig(A);
    e = diag(E);
    [e_ordered,order] = sort(e,'descend');
    V = V(:,order);
    V = V(:,2:n);
    E = diag(e_ordered(2:n));
    diff_map = E^t * V';
    diff_map = diff_map';
end