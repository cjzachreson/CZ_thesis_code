clear all
close all

As = []
Y0s = []
lambdas = []
N_on_ns = [1, 5, 10, 50, 100, 500, 1000]
n_bins =[10, 50, 100, 500, 1000, 5000, 10000]

dS_avgs = ones(numel(N_on_ns), numel(n_bins));
fit_vals = ones(numel(N_on_ns), numel(n_bins));

N_on_n_i = 0;

for N_on_n = N_on_ns
    
N_on_n_i = N_on_n_i + 1
    
    n_i = 0;
    
    for n = n_bins
        
    n_i = n_i + 1
    

        N = N_on_n * n; %number of random variables

        k = 1000; %number of distributions
        H = zeros(1, k);
        dS = zeros(1, k);

        minimum = 0;
        maximum = 1;

        %simulating sampling in time
        for t = 1:k

            vars = rand(1, N) .* (maximum - minimum) + minimum ;

            H(t) = S_dist(vars, minimum, maximum, n);

            dS(t) = (log(n) - (H(t))) ./ log(n);

        end

        dS_avg = mean(dS);
        dS_avgs(N_on_n_i, n_i) = dS_avg;

        dS_min = min(dS);

        if n_i == 1
            dS_dist_min = 0;
            dS_dist_max = max(dS);
            edges = [dS_dist_min:dS_dist_max/100:dS_dist_max];
        end


%        dS_dist = histogram(dS, edges, 'Visible', 'off');

%         figure(1)
%         plot(dS_dist.BinEdges(1:end-1), smooth(dS_dist.Values), '-')
%         drawnow
%         hold on

%         figure(2)
%         loglog(C, dS_avgs, '*')
%         drawnow


    end

%     dS_vs_n_fit = fittype('dS_vs_n_fittype(x, lambda, n_on_k, A, B)', ...
%         'problem', {'lambda', 'n_on_k', 'A', 'B'});
% 
%     ds_vs_n = fit(n_bins', dS_avgs(N_on_n_i, :)', dS_vs_n_fit,...
%         'problem', {-0.38, 1./N_on_n, 0.37, 0.045});
% 
%     coeffs =   coeffvalues(ds_vs_n);
    

    fit_vals(N_on_n_i, :) = (1/N_on_n) * (0.37 * n_bins.^(-0.38) + 0.045);

    
end

    figure(3)
     loglog(n_bins, dS_avgs', 'o')
     hold on
    loglog(n_bins, fit_vals', '--')
   

    hold off
    
    
    
    
    
    
    
    
    
    