clear all
close all

% symmetry analysis

% redimensions XY coords to polar (r, phi) and makes a histogram of density with
% respect to phi. The entropy of this distribution as a function of time
% is the symmetry metric I will use to find the time constant for evolution
% towards the isotropic state.

% requires function S_dist() to find entropy of distributions


rep_cycles = {'02', '03', '04', '05','06'...
    ,'07','08','09','10','11', '12', '13', '14'}; 
% replication cycle 

final_rep = str2double(rep_cycles{end});
first_cycle = str2double(rep_cycles{1});

dS_mat = [];

dS_max = 0;

root = ['/Users/11678505/Desktop/live/simulations/'...
    'GTR_no_stig_dt_0p001_tr_1000_GR_0.001_Pa_0.1/'];

filename = [root, 'XYTL_tr_1000_GR_0.001_Pa_0.1.txt'];

XYTL = XYTL_data(filename, 350);

output_root = [root, 'stats_vs_cycle/'];
if ~isdir(output_root)
    mkdir(output_root)
end

t_vals = [];

growth_rate = 0.001;
dl = 4;
t_rep = dl/growth_rate;
t_rec = 100;
entries_per_cycle = round(t_rep / t_rec);


for t_i = 1:numel(rep_cycles)
 
    cycle = str2double(rep_cycles{t_i})
    
    entries = ((entries_per_cycle * (cycle - 1)) : (entries_per_cycle * cycle)) + 1;
    %entries_per_cycle * cycle - 1;
    
    X_cycle = XYTL.X(:, entries);
    Y_cycle = XYTL.Y(:, entries);
    
    X_in = X_cycle(X_cycle ~= 0);
    Y_in = Y_cycle(Y_cycle ~= 0);
    
    X_cm = mean(X_in);
    Y_cm = mean(Y_in);

    rep_max = max(str2double(rep_cycles));

    max_N = 2^rep_max; 

    dims = size(XYTL.X);

    numframes = dims(2);

    phi = atan2(Y_in - Y_cm, X_in - X_cm);

    phi_0 = atan2( -Y_cm, -X_cm);

    phi = phi(phi~=phi_0) + pi;

    num_samples = numel(phi(X_in~=0));
    
    numbins = round(num_samples / 10);

    S_uni = log(numbins);
    
    dS_min = (numbins/num_samples) * ( 0.37 * numbins^(-0.38) + 0.045 );
    
    dS_tau_criteria = dS_min + dS_min * 0.05;

    [S_phi, hist_vals] = S_dist(phi, 0, 2*pi, numbins);

    dS_mat = [dS_mat, (1 - S_phi / S_uni) - dS_min];


    t_vals = [t_vals, str2double(rep_cycles{t_i})];
    
    figure(2)
    plot(t_vals, dS_mat, '-o')
    title('\DeltaS_{\phi} vs. replication cycle')
    xlabel('t_{rep}');
    ylabel('\DeltaS_{\phi}');
    hold on
    ax = gca;
    set(ax, 'FontSize', 16)
    drawnow
    
end

dS_vs_cycle_name = [output_root, 'dS_phi_vs_cycle.txt'];

dlmwrite(dS_vs_cycle_name, [dS_mat', t_vals'])
    
    



 


  












