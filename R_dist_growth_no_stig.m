clear all
close all

% expansion analysis

% redimensions XY coords to polar (r, phi) and makes a histogram of density with
% respect to r. Identify the point at which expansion becomes exponential


% ks and gs define the properties of the substratum on which the bacteria
% are moving. ks is how quickly the substratum deforms in response to
% bacterial activity, gs is the magnitude of force applied by the
% deformation gradient


growth_rate = 0.001;
dl = 4;

 rep_cycles = {'01','02','03','04','05',...
     '06', '07','08', '09','10','11','12','13','14','15'}; % replication cycle 

%  rep_cycles = {'03',...
%      '06', '09','12','15'}; % replication cycle 

%rep_cycles = {'15'};



root = ['/Users/11678505/Desktop/live/simulations/'...
    'GTR_no_stig_dt_0p001_tr_1000_GR_0.001_Pa_0.1/'];

output_root = [root, 'stats_vs_cycle/'];

if ~isdir(output_root)
    mkdir(output_root)
end


filename = [root 'XYTL_tr_1000_GR_0.001_Pa_0.1.txt'];
    
XYTL = XYTL_data(filename, 700);

dims = size(XYTL.X);

r_max = 350;
dr = 10;
r_edges = [0:dr:r_max];
A_r = pi * (r_edges + dr).^2 - pi * r_edges.^2; % area of each radial bin

A_i_avg = 5 + pi/4;


t_rep = dl/growth_rate;
t_rec = 100;

entries_per_cycle = round(t_rep / t_rec);

coverage_max = [];
FWHM = [];
r_avg = [];
t_vals = [];

legend_str = cell(1, numel(rep_cycles));
legend_hs = [];

r_cm = [];

num_frames = t_rep / t_rec;

for t_i = 1:numel(rep_cycles)
 
    cycle = str2double(rep_cycles{t_i})
    
    legend_str{t_i} = ['t_{rep} = ' rep_cycles{t_i}];
    
    entries = entries_per_cycle * cycle - 1;%((entries_per_cycle * (cycle - 1)) : (entries_per_cycle * cycle)) + 1;
  


        X_cycle = XYTL.X(:, entries);
        Y_cycle = XYTL.Y(:, entries);
  
        numframes = dims(2);

       

        X_in = X_cycle(X_cycle ~= 0);
        Y_in = Y_cycle(Y_cycle ~= 0);
        
        X_cm = mean(X_in);
        Y_cm = mean(Y_in);
        
        r_cm = [r_cm, sqrt(X_cm^2 + Y_cm^2)];
        
        X_corr = X_in - X_cm;
        Y_corr = Y_in - Y_cm;

        
        r = sqrt(X_corr.^2 + Y_corr.^2);

        r_avg = [r_avg, mean(r)];
        
        num_samples = numel(r);

        % find the number of cells in each radial bin
        figure(1)
        r_hist = histogram(r, r_edges, 'Visible', 'off');
        
        %number of cells in each bin divided by the area of the corresonding bin
        rho_hist = r_hist.Values ./ A_r(1:end-1); 
        
        coverage_hist = rho_hist .* A_i_avg;
        
        r_rad_sym_dist = [-fliplr(r_edges(2:end-1)), r_edges(1:end-1)];
        
        coverage_rad_sym_dist = ...
            [fliplr(coverage_hist(2:end)), coverage_hist];
        
%         coverage_rad_sym_dist_norm =...
%             coverage_rad_sym_dist./max(coverage_rad_sym_dist);
%         
%         gauss_norm_fittype = fittype('a1*exp(-((x-b1)/c1)^2)', 'problem', 'a1')
%         
%         gaussfit_norm = fit(r_rad_sym_dist',coverage_rad_sym_dist_norm',...
%             gauss_norm_fittype, 'problem', 1, 'StartPoint', [0, 10]) 
           
%         coverage_max = [coverage_max, max(coverage_hist)];
%         FWHM = [FWHM, 2.355 * gaussfit_norm.c1];

        

        fig = figure(2);
        h = plot(r_rad_sym_dist', coverage_rad_sym_dist', '--');
        legend_hs = [legend_hs, h];
        
        title('coverage vs r')
        xlabel('r');
        ylabel('\langle \rho \rangle');
        hold on
        %plot(gaussfit_norm, 'k-')
        lims = [0, 350, 0, 1.15];
        axis(lims)
        ax = gca;
        set(ax, 'FontSize', 16)
      
        
        if t_i == numel(rep_cycles)
           legend(legend_hs, legend_str{1:t_i})
        end
        
       drawnow
       
       cov_vs_r_name = [output_root 'cov_vs_r_cy_' rep_cycles{t_i} '.txt']
       dlmwrite(cov_vs_r_name, [coverage_rad_sym_dist', r_rad_sym_dist'])
        
       
       
        
    
    t_vals = [t_vals, str2double(rep_cycles{t_i}) * t_rep];
 
end


  
        exp_growth_fittype = fittype('a1*exp(k*x)', 'problem', 'k');
      
        k = (growth_rate / dl) * log(2) * 0.5;
     
        
         r_vs_t_fit = fit(t_vals', r_avg',...
             exp_growth_fittype, 'problem', k , 'StartPoint', [0]);
        
        exp_growth_fit_vals = 1 .* exp(k.*t_vals);
        
        D1 = 0.17;
        r1 = 1;
        ti = 1000; %time of first division event r_avg ~= 0
        r_offset = sqrt(D1 .* ti);
        r_diffusion = sqrt(D1 .* (t_vals));
        

        g = figure('Color',[1 1 1]);
        plot(t_vals, r_avg, 'ko')
        hold on
        plot(t_vals, exp_growth_fit_vals, 'r--')
        xlabel('time');
        ylabel('\langle r \rangle');
        title('r_{avg}')
        plot(t_vals - t_i, r_diffusion - r_offset, 'k--')
        set(gca, 'FontSize', 16)
        hold off
        
        
        exp_fit_name = [output_root, 'exp_fit_r_avg_vs_t.txt']; 
        diff_fit_name = [output_root, 'diff_fit_r_avg_vs_t.txt'];
        r_avg_name = [output_root, 'r_avg_vs_t.txt'];
        
        dlmwrite(exp_fit_name, [exp_growth_fit_vals', t_vals'])
        dlmwrite(diff_fit_name, [(r_diffusion - r_offset)', (t_vals - t_i)'])
        dlmwrite(r_avg_name, [r_avg', t_vals'])
        
        
        
        
        

