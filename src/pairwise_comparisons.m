%%% obtain pairwise comparisons
%   alpha :: probability of getting a pairwise comparison
%   beta  :: probability of pairwise comparison being correct
%   P     :: permutation matrix from idx --> idx_oo

function [W_oo] = pairwise_comparisons(alpha,beta,idx_oo,P)

    n = size(P,1);
    % 1 if we get this pairwise comparison, 0 otherwise
    mask_1 = logical(binornd(1,alpha,n,n));
    % 1 if pairwise comparison is incorrect, 0 otherwise
    mask_2 = logical(binornd(1,1-beta,n,n));
    % (i,j)=1 if j to the right of 1
    correct = triu(ones(n),1);
    W = correct;
    % flip incorrect entries
    W(mask_2) = ~correct(mask_2);
    % correct should be 1, incorrect -1 and diagonal 0
    W = W*2-ones(n); W(logical(eye(n))) = 0;
    % entries we do not see should be zero
    W(~mask_1) = 0;
    % enforce consistency
    W = triu(W,1);
    W = W-W';
    % permute entries to match the order of the permuted images
    W_oo = W(idx_oo,idx_oo);

end