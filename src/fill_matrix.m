function [ T_complete ] = fill_matrix( T, W )
%fill_matrix: Takes a matrix T filled with 1's and 0's 
% and uses pairwise affinities to extract information about global ordering
n = size(W,1); 

%T is the matrix of sparse pairwise comparisons, T_ij =1 if j > i -1 if i >
%j
[I,J,V ] = find(T); 
s = nnz(T); %should be equal to size of V

% make t be 1 or -1 in order to restrict attention to upper triangular part
T_complete = zeros(n);
%T complete is the matrix of -1,1s 
%Examine upper triangular part

for l=1:s
    if (V(l) == 1) %j is to right of i
        i = I(l);
        j = J(l);
        %fill row with info from W
        for k=1:n
            if k == j
                T_complete(i,k) = 1;            
            else
               if (W(i,j)> W(i,k) && W(i,k) > W(j,k)) % k i j - good
                   if T_complete(i,k) == 1
                      T_complete(i,k) = -42;
                   else
                      T_complete(i,k) = -1;
                   end
                   if T_complete(j,k) == 1
                      T_complete(j,k) = -42;
                   else
                      T_complete(j,k) = -1;
                   end
               elseif (W(i,j)> W(j,k) && W(j,k) > W(i,k)) % i j k - good
                   if T_complete(i,k) == -1
                      T_complete(i,k) = -42;
                   else
                      T_complete(i,k) = 1;
                   end
                   if T_complete(j,k) == -1
                      T_complete(j,k) = -42;
                   else
                      T_complete(j,k) = 1;
                   end
               elseif (W(j,k)> W(i,k) && W(i,k) > W(i,j)) % i k j - good
                   if T_complete(i,k) == -1
                      T_complete(i,k) = -42;
                   else
                      T_complete(i,k) = 1;
                   end
                   if T_complete(j,k) == 1
                      T_complete(j,k) = -42;
                   else
                      T_complete(j,k) = -1;
                   end
%                elseif (W(j,k)> W(i,j) && W(i,j) > W(i,k)) % i j k - good
%                    if T_complete(i,k) == -1
%                       T_complete(i,k) = -42;
%                    else
%                       T_complete(i,k) = 1;
%                    end
%                    if T_complete(j,k) == -1
%                       T_complete(j,k) = -42;
%                    else
%                       T_complete(j,k) = 1;
%                    end             
%                elseif (W(i,k)> W(i,j) && W(i,j) > W(j,k)) % k i j - good
%                    if T_complete(i,k) == 1
%                       T_complete(i,k) = -42;
%                    else
%                       T_complete(i,k) = -1;
%                    end
%                    if T_complete(j,k) == 1
%                       T_complete(j,k) = -42;
%                    else
%                       T_complete(j,k) = -1;
%                    end
%                elseif (W(i,k)> W(j,k) && W(j,k) > W(i,j)) % i k j - good
%                    if T_complete(i,k) == -1
%                       T_complete(i,k) = -42;
%                    else
%                       T_complete(i,k) = 1;
%                    end
%                    if T_complete(j,k) == 1
%                       T_complete(j,k) = -42;
%                    else
%                       T_complete(j,k) = -1;
%                    end                                                      
               end    
            end
        end
        %fill ith and jth col with info from W
        for k=1:n
            if k == i
                T_complete(k,j) = 1;            
            else             
               if (W(i,j)> W(i,k) && W(i,k) > W(j,k)) % k i j - good
                   if T_complete(k,i) == -1
                      T_complete(k,i) = -42;
                   else
                      T_complete(k,i) = 1;
                   end
                   if T_complete(k,j) == -1
                      T_complete(k,j) = -42;
                   else
                      T_complete(k,j) = 1;
                   end                    
               elseif (W(i,j)> W(j,k) && W(j,k) > W(i,k)) % i j k - good
                   if T_complete(k,i) == 1
                      T_complete(k,i) = -42;
                   else
                      T_complete(k,i) =-1;
                   end
                   if T_complete(k,j) == 1
                      T_complete(k,j) = -42;
                   else
                      T_complete(k,j) = -1;
                   end
               elseif (W(j,k)> W(i,k) && W(i,k) > W(i,j)) % i k j - good
                   if T_complete(k,i) == 1
                      T_complete(k,i) = -42;
                   else
                      T_complete(k,i) = -1;
                   end
                   if T_complete(k,j) == -1
                      T_complete(k,j) = -42;
                   else
                      T_complete(k,j) = 1;
                   end
%                elseif (W(j,k)> W(i,j) && W(i,j) > W(i,k)) % i j k - good
%                    if T_complete(k,i) == 1
%                       T_complete(k,i) = -42;
%                    else
%                       T_complete(k,i) = -1;
%                    end
%                    if T_complete(k,j) == 1
%                       T_complete(k,j) = -42;
%                    else
%                       T_complete(k,j) = -1;
%                    end              
%                elseif (W(i,k)> W(i,j) && W(i,j) > W(j,k)) % k i j - good
%                    if T_complete(k,i) == -1
%                       T_complete(k,i) = -42;
%                    else
%                       T_complete(k,i) = 1;
%                    end
%                    if T_complete(k,j) == -1
%                       T_complete(k,j) = -42;
%                    else
%                       T_complete(k,j) = 1;
%                    end
%                elseif (W(i,k)> W(j,k) && W(j,k) > W(i,j)) % 
%                    if T_complete(k,i) == 1
%                       T_complete(k,i) = -42;
%                    else
%                       T_complete(k,i) = -1;
%                    end
%                    if T_complete(k,j) == -1
%                       T_complete(k,j) = -42;
%                    else
%                       T_complete(k,j) = 1;
%                    end                                                     
               end    
            end
        end
    end   
end

T_complete(T_complete==-42) = 0;
check = T_complete+T_complete';
[ii,jj,~] = find(check);
for idx=1:length(ii)
    T_complete(ii(idx),jj(idx)) = 0;
    T_complete(jj(idx),ii(idx)) = 0;
end

end

