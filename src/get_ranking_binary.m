function [sol,order] = get_ranking_binary(W,T,lambda)
    n = size(W,1);
    
    %% DEFINE PROBLEM
    %  :: variables are t,tau(+),tau(-)
    %%% cost function
    c = [W(:); lambda*ones(2*n^2,1)];
    n_var = length(c);
    %%% linear constraints
    % t(i,j)+t(j,i)==1 :: n*(n-1)/2 constraints
    n_cons = n*(n-1)/2;
    A1 = zeros(n_cons,n_var);
    l = 1;
    for i=1:n
        for j=i+1:n
            temp = sparse([i j],[j i],[1 1],n,n);
            A1(l,:) = [temp(:)' zeros(1,n_var-n^2)];
            l = l+1;
        end
    end
    blc1 = ones(n_cons,1); buc1 = ones(n_cons,1);
    % 1<=t(i,j)+t(j,k)+t(k,i)<=2 :: n^3 constraints
%     n_cons = n*(n-1)/2*;
%     A2 = zeros(n_cons,n_var);
    n_cons = 0;
    for i=1:n
        for j=i+1:n
            n_cons = n_cons + n-j;
        end
    end
    A2 = spalloc(n_cons,n_var,3*n_cons);
    l = 1;
    for i=1:n
        for j=i+1:n
            for k=j+1:n
                A2(l,(i-1)*n+j) = 1;
                A2(l,(j-1)*n+k) = 1;
                A2(l,(k-1)*n+i) = 1;
%                 temp = sparse([i j k],[j k i],[1 1 1],n,n);
%                 A2(l,:) = [temp(:)' zeros(1,n_var-n^2)];
                l = l+1;
            end
        end
    end
    fprintf('n_cons = %d\n',size(A2,1));
    fprintf('n_cons = %d\n',n_cons);
    blc2 = ones(size(A2,1),1); buc2 = 2*ones(size(A2,1),1);
    % tau(+)_ij-tau(-)_ij-t_ij=-t(hat)_ij :: sames as number of edges
    l = 1;
    [ii,jj,s] = find(T);
    n_cons = length(ii);
    A3 = zeros(n_cons,n_var);
    for idx=1:n_cons
        i = ii(idx); j = jj(idx);
        temp = [sparse(i,j,-1,n,n) sparse(i,j,1,n,n) sparse(i,j,-1,n,n)];
        A3(l,:) = temp(:)';
        l = l+1;
    end
    blc3 = -s; buc3 = -s;
    % all together
    A = sparse([A1;A2;A3]);
    blc = [blc1;blc2;blc3];
    buc = [buc1;buc2;buc3];
    %%% upper and lower bounds
    blx = zeros(n_var,1); bux = ones(n_var,1);
    
    %% SOLVE
    res = msklpopt(c,A,blc,buc,blx,bux);
    
    %% FIND RANKING ORDER
    sol = reshape(res.sol.itr.xx(1:n^2),n,n);
    gt = sol>=0.5;  % 1 if j>i
    csums = sum(gt,1);
    [~,order] = sort(csums);
    order = order';
    
end