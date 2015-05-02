clear; clc;
%%% parameters
alpha   = 0.01;       % fraction of pairwise comparisons
beta    = 0.90;       % probability of correct pairwise comparisons

%%% get images
range = 120;
imgs = image_reader('zebrafish',range);
%%% convert the uint8 pixels to doubles
imgs = double(imgs);

%%% number of images, should equal range
n = size(imgs,3);

%%% permute images
P = eye(n); P = P(randperm(n),:);
idx = 1:n;
idx_oo = idx * P;       % takes 1:n to idx_oo
idx_rv = idx * P';      % takes idx_oo to 1:n
imgs_oo = imgs(:,:,idx_oo);

%%% calculate image affinities (weights)
A = gaussian_kernel_weights(imgs_oo);
%%% matrix of degrees: (completely connected)
D = diag(sum(A));
L = inv(D)*A; 

%%% get pairwise comparison
W = pairwise_comparisons(alpha,beta,idx_oo,P);

%%% 1. (vector) diffusion maps
[V,E] = eig(L);
[E order] = sort(diag(E),'descend');  % sort eigenvalues in descending order
V = V(:,order);
%Discard the first eiegenvector of all one's
V = V(:,2:n);
E= E(2:n);

t=2; %how do we define this?

%R1 embeddding
x_pos = E(1).^t*sqrt(inv(D(2,2)))*V(:,1);
figure;
plot(x_pos,zeros(1,n),'-o');
%R2 embedding
y_pos = E(2).^t*sqrt(inv(D(3,3)))*V(:,2);
figure;
plot(x_pos,y_pos,'-o')

%%% 2. ranking + pairwise comparisons


%%% 3. ranking + pairwise comparisons + time stamps 

%%% 4. doubly stochastic relaxation

%%% 5. local nonconvex relaxation

%%% 6. spectral relaxation

%%% 7. esdp relaxation