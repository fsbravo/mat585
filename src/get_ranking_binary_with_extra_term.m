function [sol,order] = get_ranking_binary_with_extra_term(T,E)
    n = size(T,1);
    T(T<0) = 0;
    
    %% PARAMETERS
    lambda = 1;
    
    %% DEFINE PROBLEM
    %  :: variables are t,tau(+),tau(-),tau3(+),tau3(-)
    %%% cost function
    [ii,jj,ss] = find(triu(T,1));
    [iiE,jjE,ssE] = find(triu(E,1)); n_cons_3 = n*length(ssE);
    c = [zeros(n^2,1); ones(2*length(ss),1); lambda*ones(2*n_cons_3,1)];
    n_var = length(c);
    %%% linear constraints
    % t(i,j)+t(j,i)==1 :: n*(n-1)/2 constraints
    n_cons = n*(n-1)/2;
    A1 = spalloc(n_cons,n_var,2*n_cons);
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
    n_cons = length(ii);
    A3 = spalloc(n_cons,n_var,3*n_cons);
    for idx=1:n_cons
        i = ii(idx); j = jj(idx);
        temp = sparse(i,j,-1,n,n);
        A3(l,:) = [temp(:)',sparse(1,idx,1,1,length(ss)),sparse(1,idx,-1,1,length(ss)+2*n_cons_3)];
        l = l+1;
    end
    blc3 = -ss; buc3 = -ss;
    % tau3(+)_ijk-tau(-)_ijk-t_ik+t_ij=0;
    A4 = spalloc(n_cons_3,n_var,4*n_cons_3);
    size(A4)
    l = 1;
    for idx=1:length(ssE)
        i = iiE(idx); j = jjE(idx);
        for k=1:n
            v_idx = (idx-1)*n+k;
            temp = sparse([i,j],[k,k],[-1,+1],n,n);
            A4(l,:) = [temp(:)',zeros(1,2*length(ss)),...
                sparse(1,v_idx,1,1,n_cons_3),sparse(1,v_idx,1,1,n_cons_3)];
            l = l+1;
        end
    end
    blc4 = zeros(n_cons_3,1); buc4 = zeros(n_cons_3,1);
    % all together
    A = sparse([A1;A2;A3;A4]);
    blc = [blc1;blc2;blc3;blc4];
    buc = [buc1;buc2;buc3;blc4];
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