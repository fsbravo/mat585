function [t,D] = get_ranking_base(W,T)

n = size(W,1);

cvx_begin
    cvx_solver mosek
    cvx_precision low
    variable t(n,1)
    variable D(n,n)
    minimize( ranking_obj(t,D,T,W) );
    subject to
        D(:)>=1
cvx_end

end