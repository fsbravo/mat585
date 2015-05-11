%%% augmented diffusion maps with time stamps
%   W :: adjacency matrix
%   T :: matrix of pairwise comparisons
%   t_hat :: vector of noisy time stamps
%   lambda :: weight of pairwise term
%   gamma :: weight of time stamps
function [t,d] = get_ranking_base_time(W,T,t_hat,lambda,gamma)
D = diag(sum(W,2));
L = D-W;
n = size(W,1);

cvx_begin
    cvx_solver mosek
    cvx_precision low
    variable t(n,1)
    variable d(n,n)
    minimize( ranking_obj_time(t,d,T,L,t_hat,lambda,gamma) );
    subject to
        d(:)>=1;
        t'*ones(n,1) == 0;
%         t(n)-t(1) == n;
cvx_end

end