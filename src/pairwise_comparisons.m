%%% obtain pairwise comparisons
%   alpha :: probability of getting a pairwise comparison
%   beta  :: probability of pairwise comparison being correct
%   P     :: permutation matrix from idx --> idx_oo

function [W] = pairwise_comparisons(alpha,beta,P)

    n = size(P,1);
    W = zeros(n);
    mask = binornd(1,alpha,n,n);
    correct = triu(ones(n),-1) + triul(ones(n),-1);

end