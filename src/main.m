clear; clc;
%%% parameters
alpha   = 1.00;       % fraction of pairwise comparisons
beta    = 1.00;       % probability of correct pairwise comparisons

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
W = gaussian_kernel_weights(imgs);
%%% matrix of degrees: (completely connected)
D = diag(sum(W));

%%% get pairwise comparison
T = pairwise_comparisons(alpha,beta,idx,P);

%%% 1. (vector) diffusion maps
t=1; %how do we define this?
diff_map = full_diffusion_map(W,t);
x_pos = diff_map(:,1);
y_pos = diff_map(:,2);
%R1 embeddding
figure;
plot(x_pos,zeros(1,n),'-o');
%R2 embedding
figure;
plot(x_pos,y_pos,'-o')

%%% 2. ranking + pairwise comparisons
order_2 = get_ranking_base(W,T);

%%% 3. ranking + pairwise comparisons + time stamps 

%%% 4. doubly stochastic relaxation

%%% 5. local nonconvex relaxation

%%% 6. spectral relaxation

%%% 7. esdp relaxation