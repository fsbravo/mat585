function [val] = ranking_obj(t,D,T,W)
    n = length(t);
    val = 0;
    lambda = 1;
    for i=1:n
        for j=1:n
            val = val + abs(t(i)-t(j))*W(i,j) + lambda*abs(t(j)-t(i)-D(i,j)*T(i,j));
        end
    end
end