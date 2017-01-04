% Script used to calculate data in figure 4.22(f) lambda vs gamma

close all
clear all



gs = {'0.001','0.25','0.50','0.60','0.70','0.80', '0.90' ,'1.0','1.1','1.2','1.3','1.4','1.5','1.6'};

thresholds = {'8', '10', '16', '20', '32'} ;


output_mat = NaN(numel(gs), numel(thresholds));


listindex = 0;

scaling_coeffs = [];
confidence_interval = [];

for gs_i = 1:length(gs)
    
fluct_n = [];
equilibrium_n = [];

for th_i = 1:length(thresholds);
for gl_i = 1:length(gl)
for s_i = 1:length(state)
for tr_i = 1:length(Trev)

listindex = listindex + 1;
    
% root = ['/Users/11678505/Desktop/PA_sim_data/PB_followup/tf_100k/'];
root = ['/Users/11678505/Desktop/PA_sim_data/PB_gl_0_all/'];
XYTL_filename = [root ...
    '/XYTL/XYTL_T_rev_' Trev{tr_i} '_test_g_s_' gs{gs_i} '.txt'];

L = 160

l = L / str2double(thresholds{th_i});

output_path = [root 'density_fluct/']
if ~isdir(output_path)
    mkdir(output_path)
end

L = 160

XYTL = XYTL_data(XYTL_filename, L/2);

X = XYTL.X;
Y = XYTL.Y;

dims = size(XYTL.l);

clear XYTL

numframes = dims(2) - 1;

N = dims(1);

IDs = 1:N;

grid_res = l;


 XI_max = L / grid_res;
 YI_max = L / grid_res;
 
ti = numframes - numframes/5;

loc_dens = zeros(XI_max*YI_max, numel(ti:numframes));
%exchange_grid = zeros(N, numframes);

 XI = floor(X / grid_res) + (L / grid_res / 2) + 1;
 YI = abs(ceil(Y / grid_res) - (L / grid_res / 2)) + 1;
% 
 XI(X == L/2) = XI_max;
 YI(Y == - L/2) = YI_max;
%fill structure grid
for t = ti:numframes
    t
    for xi = 1:XI_max
        for yi = 1:YI_max
            lind = (xi - 1) * YI_max + yi;
            loc_dens( lind, t ) = numel(IDs(XI(:, t) == xi & YI(:, t) == yi));
        end
    end
            
    
    
end

frames = ti:numframes-1;

loc_dens_stst = loc_dens(:,frames);

avg_n = N/L^2 * l^2; %what we expect for equilibrium 

fluct_n_l = sqrt(mean(mean((loc_dens_stst - avg_n).^2)));

fluct_n = [fluct_n, fluct_n_l];

equilibrium_n = [equilibrium_n, avg_n];


output_mat(gs_i, th_i) = fluct_n_l;

   

   %dlmwrite([output_path 'r_' thresh_txt '_local_density_gs_' gs{gs_i} '.txt'], h') 

end
end
end
end


lin_fit = fittype('lin_fit(x, m, b)');

figure(4)

plot(log(equilibrium_n), log(fluct_n), '-*')

hold on

 fit_f1 = fit(log(equilibrium_n'), log(fluct_n'), lin_fit, 'StartPoint', [0, 0.5] );
% 
% 
 b1 = coeffvalues(fit_f1);
 ci = confint(fit_f1, 0.95); fitvals =   exp(b1(1)) * equilibrium_n .^ b1(2);
 ci_lambda = (ci(2, 2) - ci(1, 2)) / 2;
 scaling_coeffs = [scaling_coeffs, b1(2)]
 confidence_interval = [confidence_interval, ci_lambda]
 
 figure(2) 
 loglog(equilibrium_n, fluct_n, 'o')
 hold on
 loglog(equilibrium_n, fitvals, 'k--')
 
drawnow


end

 gs_axis = [0.001, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6];
% 
% 
% figure(3)
% plot(gs_axis, scaling_coeffs, '*')
 dlmwrite([output_path 'scaling_coefficients.txt'], [scaling_coeffs', gs_axis'])
 dlmwrite([output_path 'confidence_intervals.txt'], confidence_interval')

 dlmwrite([output_path  'y_gs_x_l_z_delta_n.txt'], [output_mat])








