clear all
close all

% symmetry analysis

% finds avrage anisotropy in local neighborhoods of increasing radius


rep_cycles = {'11'}%{'02', '03', '04', '05','06'...
    %,'07','08','09','10','11'}%, '12', '13', '14'}; 
% replication cycle 


scales = [3:1:20, 21:5:50, 51:10:100, 101:50:300];


gs = {'0.5'};
ks = {'0.05', '0.001'};


root = ['/Users/11678505/Desktop/live/simulations/Propagating_anisotropy/'...
    'furrowing_only/PA_P_0.5_L100_dt_0p005/' ]


dS_max = 0;

for g_i = 1:numel(gs)
for k_i = 1:numel(ks)

output_root = [root, 'dS_vs_cycle/'];
if ~isdir(output_root)
    mkdir(output_root)
end

r_max = 350;


growth_rate = 0.001;
dl = 4;
t_rep = dl/growth_rate;
t_rec = 100;
entries_per_cycle = round(t_rep / t_rec);


for t_i = 1:numel(rep_cycles)
  
    
    XYTL_root = [root 'PA_P_0.5_L100_dt_0p005_ks_' ks{k_i} '_gs_' gs{g_i} '/'];
    
    XYTL_filename = [XYTL_root,...
        'XYTL_trep_' rep_cycles{t_i} '_ks_' ks{k_i} '_gs_' gs{g_i} '.txt'];

    

    XYTL = XYTL_data(XYTL_filename, 350);
 
    X_cycle = XYTL.X;
    Y_cycle = XYTL.Y;
    
    X_in = X_cycle(X_cycle ~= 0);
    Y_in = Y_cycle(Y_cycle ~= 0);
    
    N = numel(X_in);
    
    IDs = 1:N;

    num_i = 1000;%round(N / 1000);
    
    sampled_cells = [];
    
    n_unique = numel(unique(sampled_cells));
    
    while n_unique < num_i
    
    sampled_cells = [unique(sampled_cells), ceil(rand(1, num_i - n_unique) * N)];
    
    n_unique = numel(unique(sampled_cells));
    end
    
    
    dS_mat = zeros(1, numel(scales));
    
    for l_i = 1:numel(scales)
    
    l = scales(l_i)

    dS_loc = zeros(1, num_i);
    
    for n_i = 1:num_i

    i = sampled_cells(n_i);
    
    X_i = X_in(i);
    Y_i = Y_in(i);

    dXij = X_in - X_i;
    dYij = Y_in - Y_i;
    
    dij = sqrt(dXij.^2 + dYij.^2);
    
    
    hood = IDs(dij < l);
    hood = hood(hood ~= i);
    
%     figure(4)
%     plot(X_in(hood), Y_in(hood), '.')
%     hold on
%     plot(X_i, Y_i, 'o')
%     hold off
%     drawnow
    
    phi = atan2(Y_in(hood) - Y_i, X_in(hood) - X_i) + pi;

    num_samples = numel(hood);
    
    numbins = round(num_samples / 10);

    S_uni = log(numbins);
    
    dS_min = (numbins/num_samples) * ( 0.37 * numbins^(-0.38) + 0.045 );

    [S_phi, hist_vals] = S_dist(phi, 0, 2*pi, numbins);

    dS_loc(n_i) = (1 - S_phi / S_uni) - dS_min;
    
    
    end
    
    dS_mat(l_i) = mean(dS_loc);
    
    end
    
    figure(2)
    semilogx(scales, dS_mat, 'o')
    title('\DeltaS_{\phi} vs. replication cycle')
    xlabel('r');
    ylabel('\DeltaS_{\phi}');
    hold on
    ax = gca;
    set(ax, 'FontSize', 16)
    drawnow
   
    dlmwrite([output_root 'cycle_' rep_cycles{t_i} '_ks_' ks{k_i} '_gs_'...
        gs{g_i} '.txt'], [dS_mat', scales'])
    
   
end
    
    

end
end

 


  












