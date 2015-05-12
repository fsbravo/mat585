clear; 
clc;
%Turn plots on/off
PLOT_SWITCH = 0;

%For testing a range of sampling ratios
NALPHA = 5;
NBETA = 3;
%%% Set actual parameters/ranges here as vectors
alpha   = [0.04 0.08 0.12 0.16 0.20];       % fraction of pairwise comparisons
beta    = [0.95 0.98 1.00];       % probability of correct pairwise comparisons

%Pairwise comparisons are randomly sampled. 
%Should run multiple random experiments and average results 
NEXPERIMENTS = 10;

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

%%% create matrix of errros
metrics = zeros(5,4,NEXPERIMENTS,NALPHA,NBETA); %number of solution methods \times number of metrics

%%% matrix 
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

%1. Diffusion map sort quality
t=1; %how do we define this?
diff_map = full_diffusion_map(W,t);
x_pos = diff_map(:,1);
y_pos = diff_map(:,2);
[x_sorted, ord_x] = sort(x_pos);
if (PLOT_SWITCH)
    h = figure(3);
    set(h,'name','Diffusion map ranking','numbertitle','off');
    scatter(1:range,ord_x);
    title('Quality of ranking');
    xlabel('True image rank');
    ylabel('Diffusion map rank');
end
metrics(1,:,1,1,1) = rank_metrics(ord_x,n);
for i=2:NEXPERIMENTS
    for j=1:NALPHA
        for k=1:NBETA
            metrics(1,:,i,j,k) = metrics(1,:,1,1,1);
        end
    end
end
% %R1 embeddding
% figure(1);
% plot(x_pos,zeros(1,n),'-o');
% title('Diffusion map embedding in \bf{R}^1','Interpreter','tex');
% %R2 embedding
% figure(2);
% plot(x_pos,y_pos,'-o')
% title('Diffusion map embedding in \bf{R}^2','Interpreter','tex');

%%%%%%%%%-------------BEGIN LOOP OF PAIRWISE COMPARISON SAMPLES------------
for i=1:NALPHA
    for j=1:NBETA
        for k=1:NEXPERIMENTS   
            %%% get new pairwise comparison
            T = pairwise_comparisons(alpha(i),beta(j),idx,P);
            T_complete = fill_matrix(T,W);

            %%% 2. ranking + pairwise comparisons
            [t2,D2] = get_ranking_base(W,T,0.01);
            [~, ord_t2] = sort(t2); 
            metrics(2,:,k,i,j) = rank_metrics(ord_t2,n);
            if (PLOT_SWITCH)
                h = figure(4);
                set(h,'name','Ranking + pairwise comparisons','numbertitle','off');
                %these are multicolored
                scatter(1:range,ord_t2);
                title('Quality of ranking');
                ylabel('dm+pairwise map rank');
                xlabel('True image rank');
            end

            %%% 3. ranking + pairwise comparisons + time stamps
            t_hat = 1:range;
            t_hat = t_hat' - mean(t_hat);
            sigma = 1;
            t_hat = t_hat + sigma*normrnd(0,1,range,1);
            gamma = 1;
            [t3,D3] = get_ranking_base_time(W,T,t_hat,0.01,1);
            [~, ord_t3] = sort(t3);
            metrics(3,:,k,i,j) = rank_metrics(ord_t3,n);
            if (PLOT_SWITCH)
                h = figure(5);
                set(h,'name',' Ranking + pairwise compoarisons + time stamps','numbertitle','off');
                %these are multicolored
                scatter(1:range,ord_t3);
                title('Quality of ranking');
                ylabel('dm+pairwise+time rank');
                xlabel('True image rank');
            end

            %%% 4. binary ranking method - formulate as linear program
            Tb = T; Tb(Tb<0) = 0; 
            [res,ord_b] = get_ranking_binary(Tb);
            metrics(4,:,k,i,j) = rank_metrics(ord_b,n);
            if (PLOT_SWITCH)
                h = figure(6);
                set(h,'name','Binary ranking method','numbertitle','off');
                scatter(1:range,ord_b);
                title('Quality of ranking');
                ylabel('binary rank');
                xlabel('True image rank');
            end
            
            %%% 5. binary ranking method - formulate as linear program
            Tb = T_complete; Tb(Tb<0) = 0;
            [res,ord_b] = get_ranking_binary(T_complete);
            metrics(5,:,k,i,j) = rank_metrics(ord_b,n);
            if (PLOT_SWITCH)
                h = figure(6);
                set(h,'name','Binary ranking method','numbertitle','off');
                scatter(1:range,ord_b);
                title('Quality of ranking');
                ylabel('binary rank');
                xlabel('True image rank');
            end
        end
        save temp
    end
end