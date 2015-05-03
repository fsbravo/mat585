clear; clc;
%%% parameters
alpha   = 0.50;       % fraction of pairwise comparisons
beta    = 1.00;       % probability of correct pairwise comparisons

%%% get images
range = 120;
imgs = image_reader('zebrafish',range);
%%% convert the uint8 pixels to doubles
imgs = double(imgs);
imgs_mat = zeros(size(imgs,1)*size(imgs,2),range);
for i=1:range
    img = imgs(:,:,i);
    imgs_mat(:,i) = img(:);
end
clear img;

%%% number of images, should equal range
n = size(imgs,3);

%%% permute images
P = eye(n); P = P(randperm(n),:);
idx = 1:n;
idx_oo = idx * P;       % takes 1:n to idx_oo
idx_rv = idx * P';      % takes idx_oo to 1:n
imgs_oo = imgs(:,:,idx_oo);
imgs_mat_oo = imgs(:,idx_oo);

%%% calculate image affinities (weights)
[W,distances]= gaussian_kernel_weights(imgs_mat,0.25);
%%% matrix of degrees
D = diag(sum(W));

%%% get pairwise comparison
T = pairwise_comparisons(alpha,beta,idx,P);

%%% 1. (vector) diffusion maps
t=1; %how do we define this?
diff_map = full_diffusion_map(W,t);
x_pos = diff_map(:,1);
y_pos = diff_map(:,2);
%R1 embeddding
figure(1);
plot(x_pos,zeros(1,n),'-o');
title('Diffusion map embedding in \bf{R}^1','Interpreter','tex');
%R2 embedding
figure(2);
plot(x_pos,y_pos,'-o')
title('Diffusion map embedding in \bf{R}^2','Interpreter','tex');

%Check sort quality
figure(3);
[x_sorted, ord] = sort(x_pos);
scatter(ord,1:120);
title('Quality of ranking');
xlabel('Diffusion map rank');
ylabel('True image rank');


%%% 2. ranking + pairwise comparisons
n_imgs = 10;
start = 1;
idx = start:start+n_imgs-1;
[t,D] = get_ranking_base(W(idx,idx),T(idx,idx))

%%% 3. ranking + pairwise comparisons + time stamps 

%%% 4. doubly stochastic relaxation

%%% 5. local nonconvex relaxation

%%% 6. spectral relaxation

%%% 7. esdp relaxation