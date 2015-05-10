function [t,d] = get_ranking_base(W,T,lambda)
D = diag(sum(W,2));
L = D-W;
n = size(W,1);

cvx_begin
    cvx_solver mosek
    cvx_precision low
    variable t(n,1)
    variable d(n,n)
    minimize( ranking_obj(t,d,T,L,lambda) );
    subject to
        d(:)>=1;
        t'*ones(n,1) == 0;
        t(n)-t(1) == n;
cvx_end

end