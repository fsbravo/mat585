clear; clc;
%%% parameters
alpha   = 0.01;       % fraction of pairwise comparisons
beta    = 0.90;       % probability of correct pairwise comparisons

%%% get images
% call zebrafish_reader.m
n = size(imgs,3);

%%% permute images
P = eye(n); P = P(randperm(n),:);
idx = 1:n;
idx_oo = idx * P;
imgs_oo = imgs(:,:,idx_oo);

%%% calculate image affinities
A = gaussian_kernel_weights(imgs_oo);

%%% get pairwise comparison
W = pairwise_comparisons(alpha,beta,P);

%%% 1. (vector) diffusion maps

%%% 2. ranking + pairwise comparisons

%%% 3. ranking + pairwise comparisons + time stamps 

%%% 4. doubly stochastic relaxation

%%% 5. local nonconvex relaxation

%%% 6. spectral relaxation

%%% 7. esdp relaxation