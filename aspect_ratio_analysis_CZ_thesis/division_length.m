% finds cell division events and plots histograms of length at division and
% just after division

clear all
close all

root = ['/Users/11678505/Desktop/live/Experimental_data/WT_behind_rafts/'];

L = dlmread([root, 'MajorAxisLength/MajorAxisLength.txt']);
W = dlmread([root, 'MinorAxisLength/MinorAxisLength.txt']);
X = dlmread([root, 'Centroid/X.txt']);
Y = - dlmread([root, 'Centroid/Y.txt']);
theta = dlmread([root, 'Orientation/Orientation.txt']) - 90;



IN_mat_L = sign(L);
IN_mat_XY = ones(size(X)) ;

AR_initial = [];
W_avg = [];

L_avg = mean(L(~isnan(L)));

Xmin = min(min(X));
X = X - Xmin;
Xmin = min(min(X));
Xmax = max(max(X));


Ymin = min(min(Y));
Y = Y - Ymin;
Ymin = min(min(Y));
Ymax = max(max(Y));

dims = size(L);

N = dims(1);
tf = dims(2);

IDs = 1:N;
IDs_in = IDs;


theta_in = theta;
X_in = X;
Y_in = Y; 
L_in = L;
W_in = W;

entry_point = zeros(1, numel(IDs_in));
exit_point = zeros(1, numel(IDs_in));

for i = IDs
    L_i = L(i, :);
    entry_point(i) = min(find(~isnan(L_i)));
    exit_point(i) = max(find(~isnan(L_i)));
end

L_pre_div = []
L_post_div = []

W_avg = []

for i = IDs_in

    L_i = L_in(i, :);
    
    if sum(~isnan(L_i)) > 100
    
    W_i = W_in(i, :);
    X_i = X_in(i, :);
    Y_i = Y_in(i, :);
    theta_i = theta_in(i, :);

        if L_i(exit_point(i)) < (L_i(exit_point(i) - 1) * 0.7) &&...
                L_i(exit_point(i)-1) > 50

        L_pre_div =  [L_pre_div, L_i(exit_point(i) - 1)]; 
        L_post_div = [L_post_div, L_i(exit_point(i))];
        W_avg = [W_avg, mean(W_i(~isnan(W_i)))];
        end
    
    end
end

edges = 1.5:0.25:8;

h1 = histogram(L_pre_div./W_avg, edges);
hold on
h2 = histogram(L_post_div./W_avg, edges);
hold off

AR_post_div = h2.Values;
AR_pre_div = h1.Values;

x_axis = edges(1:end-1) + 0.125;

dlmwrite([root 'AR_pre_post_div/AR_post_div_hist.txt'], [AR_post_div', x_axis'])
dlmwrite([root 'AR_pre_post_div/AR_pre_div_hist.txt'], [AR_pre_div', x_axis'])

dlmwrite([root 'AR_pre_post_div/L_post_div.txt'], L_post_div')
dlmwrite([root 'AR_pre_post_div/L_pre_div.txt'], L_pre_div')
dlmwrite([root 'AR_pre_post_div/W_avg.txt'], W_avg')


