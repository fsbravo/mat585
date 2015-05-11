function [val] = ranking_obj_time(t,d,T,L,t_hat,lambda,gamma)
    n = length(t);
    val = t'*L*t;
    [ii,jj,~] = find(T);
    for idx=1:length(ii)
        i = ii(idx);
        j = jj(idx);
        val = val + lambda*abs(t(j)-t(i)-d(i,j)*T(i,j));
    end
    val = val + gamma*sum(abs(t-t_hat));
end