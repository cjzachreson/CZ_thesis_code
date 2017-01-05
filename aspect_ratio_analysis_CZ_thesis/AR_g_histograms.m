clear all
close all

% modify root to the location of your data

type = 'WT_behind_rafts'

root = ['/Users/11678505/Desktop/live/Experimental_data/' type '/'];

L = dlmread([root, 'MajorAxisLength/MajorAxisLength.txt']);
W = dlmread([root, 'MinorAxisLength/MinorAxisLength.txt']);
X = dlmread([root, 'Centroid/X.txt']);
Y = - dlmread([root, 'Centroid/Y.txt']);
theta = dlmread([root, 'Orientation/Orientation.txt']);

IN_mat_L = sign(L);

L_cor = [];
W_cor = [];
X_cor = [];
Y_cor = [];
theta_cor = [];

growth_rates = [];
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

IDs_in = [];

t_min = 200;
break_threshold = 5;

% number of frames over which a cell must be tracked in order for its ID
% to be included in the list of IDs_in

for i = 1:N
    frames_out = numel(find(isnan(L(i, :))));
    if frames_out < (tf - t_min)    
        IDs_in = [IDs_in, i];
    end
end


for i = IDs_in
    
    X_i = X(i, :);
    Y_i = Y(i, :);
    edge_frames = find(X_i < L_avg | X_i > (Xmax - L_avg) | ...
                       Y_i < L_avg | Y_i > (Ymax - L_avg));
                   
    IN_mat_XY(i, edge_frames) = NaN;
    IN_mat_L(i, edge_frames) = 0;
                   
 end

% X = X.*IN_mat;
% Y = Y.*IN_mat;
L = L.*IN_mat_L;
W = W.*IN_mat_L;

dlmwrite('IDs_in.txt', IDs_in);

for i = IDs_in

    L_i = L(i, :);
    W_i = W(i, :);
        
    
    dL_i = [0, 0, diff(L_i), 0];
    
    L_i_pad = [0, L_i, 0];
    W_i_pad = [0, W_i, 0];
    
    
    breaks = [1, find(abs(dL_i) > break_threshold), tf+1];

    n_breaks = numel(breaks);
    
    
    stop = true;
    
    
   for j = 1:(n_breaks - 1)
  
        pre_break = 1:breaks(j);
        
        if j ~= n_breaks
            
        between_breaks = (breaks(j)+1):(breaks(j+1)-1);
        
        else
            between_breaks = (breaks(j)+1):(tf+1);
        end
        
        post_break = (breaks(j+1)+1):(tf + 1);
        
        
        L_i_break = ...
            [NaN .* L_i_pad(pre_break), ...
            L_i_pad(between_breaks), ...
            NaN .* L_i_pad(post_break)];
        
        W_i_break = ...
            [NaN .* W_i_pad(pre_break), ...
            W_i_pad(between_breaks), ...
            NaN .* W_i_pad(post_break)] ;     
        

        

        if sum(~isnan(L_i_break)) > t_min ...
                && sum(L_i_break(~isnan(L_i_break))) ~= 0 ...
                && stop
        
%                   figure(1)
%               plot(AR_break(2:end-1))
%               drawnow
            

            AR_break = (L_i_break) ./...
            (mean(W_i_break(~isnan(W_i_break))));

            L_cor = [L_cor; L_i_break(2:end-1)];
            W_cor = [W_cor; W_i_break(2:end-1)];
            W_avg = [W_avg; mean(W_i_break(end))];
            
            AR_tst = AR_break(~isnan(AR_break));
            AR_tst = AR_tst(2:end-1);


            fit_AR = polyfit([1:numel(AR_tst)]', AR_tst', 1);

            fit_vals = fit_AR(1).*[1:numel(AR_tst)]+fit_AR(2);


            growth_rate = (fit_AR(1)/2);

            growth_rates = [growth_rates; growth_rate]; %in um/s 
            AR_initial = [AR_initial, fit_vals(1)];
            stop = false;
            
        
        end
        

   end
   

   
   
end



%         figure(5)
%         h = histogram(growth_rates, 'Normalization', 'probability');
%         x = h.BinEdges;
%         dx = h.BinWidth;
%         x = x + dx/2;
%         dlmwrite([type '_GR_hist_vals.txt'], [h.Values', x(1:end-1)'])
%         dlmwrite([type '_GR.txt'], growth_rates)
%         
%         figure(6)
%         scatter(AR_initial, growth_rates)
        
        
        figure(7)
        AR_edges = 2:0.5:8;
        hh = histogram(AR_initial, AR_edges, 'Normalization', 'probability') ;      
        x = hh.BinEdges;
        dx = hh.BinWidth;
        x = x + dx/2;
        dlmwrite([type '_AR_hist_vals.txt'], [hh.Values', x(1:end-1)'])
        dlmwrite([type '_AR.txt'], AR_initial')
        
%         figure(8)
%         hhh = histogram(W_avg,'Normalization', 'probability');
%         x = hhh.BinEdges;
%         dx = hhh.BinWidth;
%         x = x + dx/2;
%         dlmwrite([type '_W_hist_vals.txt'], [hhh.Values', x(1:end-1)'])
%         dlmwrite([type '_W_Avg.txt'], W_avg)
      
        
        average_aspect_ratio = mean(AR_initial)
        std_aspect_ratio = std(AR_initial)
        average_growth_rate = mean(growth_rates)
        std_growth_rate = std(growth_rate)
        
        numel(AR_initial)
