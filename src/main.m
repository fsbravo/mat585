clear; clc;
%%% parameters
alpha   = 0.10;       % fraction of pairwise comparisons
beta    = 1.00;       % probability of correct pairwise comparisons

%%% get images
range = 120;
[imgs, nchannels] = image_reader('zebrafish',range);
%%% convert the uint8 pixels to doubles
imgs = double(imgs);
imgs_mat = zeros(size(imgs,1)*size(imgs,2)*nchannels,range);
for i=1:range
    img = imgs(:,:,:,i);
    imgs_mat(:,i) = img(:); %vectorize images. not sure if this is 
    % appropriate in RGB 3 channel case
end
clear img;

%%% number of images, should equal range
n = range; 

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
% %R1 embeddding
% figure(1);
% plot(x_pos,zeros(1,n),'-o');
% title('Diffusion map embedding in \bf{R}^1','Interpreter','tex');
% %R2 embedding
% figure(2);
% plot(x_pos,y_pos,'-o')
% title('Diffusion map embedding in \bf{R}^2','Interpreter','tex');

% check sort quality
figure(3);
[x_sorted, ord_x] = sort(x_pos);
scatter(1:range,ord_x);
title('Quality of ranking');
xlabel('Diffusion map rank');
ylabel('True image rank');
sum(abs(ord_x-[1:range]'))

%%% 2. ranking + pairwise comparisons
[t,D] = get_ranking_base(W,T,0.01);
[t_sorted, ord_t] = sort(t); hold on;
scatter(1:range,ord_t);
title('Quality of ranking');
xlabel('Diffusion map rank');
ylabel('True image rank');
sum(abs(ord_t-[1:range]'))

%%% 3. ranking + pairwise comparisons + time stamps 

%%% 4. binary ranking method - formulate as linear program
%Tb = T; Tb(Tb<0) = 0; 
%[res,A,blc,buc] = get_ranking_binary(W,Tb,1);

%%% 5. doubly stochastic relaxation

%%% 6. local nonconvex relaxation

%%% 7. spectral relaxation

%%% 8. esdp relaxation