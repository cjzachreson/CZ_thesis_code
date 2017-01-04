%script used to determine density fluctuations in Fig. 4.6 (collective behavior of twitching rods)

close all
clear all

P = {'0.1', '0.3', '0.6', '1'}
Trev = {'1000'};
cov = {'0.05', '0.1', '0.15', '0.2', '0.25', '0.3'};
intervals = {'50000'}%,'20000','30000','40000','50000'};

thresholds = {'24'} ;

listindex = 0;

scaling_coeffs = [];
confidence_interval = [];

for int_i = 1:length(intervals)
    
  legend_str = cell(1, numel(P)+1);
  legend_handles = [];
    
for P_i = 1:length(P)
    
 legend_str{P_i} = ['P_{s} = ' P{P_i}];
    
 d_rho_vs_Cov_P = [] 
 d_rho_eq_vs_Cov_P = [];
    
for cov_i = 1:length(cov)
    
  
fluct_n = [];
equilibrium_n = [];
d_rho = [];
eq_d_rho = [];

for th_i = 1:length(thresholds);

listindex = listindex + 1;
    

root = ['/Users/11678505/Desktop/live/simulations/TR_PB_no_stig/TR_polar/'];
XYTL_filename = [root, 'TR_polar_cov_' cov{cov_i} '_P_' P{P_i} ...
    '/XYTL_trep_' intervals{int_i} '_cov_' cov{cov_i} '_P_' P{P_i} '.txt'];

L = 240

l = L / str2double(thresholds{th_i});

output_path = [root 'density_fluct/']
if ~isdir(output_path)
    mkdir(output_path)
end

L = 240

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
 
ti = 1;

loc_N_interval = zeros(XI_max*YI_max, numel(ti:numframes));
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
            loc_N_interval( lind, t ) = numel(IDs(XI(:, t) == xi & YI(:, t) == yi));
        end
    end
             
end

frames = ti:numframes-1;

loc_N_stst = loc_N_interval(:,frames);

loc_avg_n = N/L^2 * l^2; %what we expect for equilibrium 

fluct_n = [fluct_n, sqrt(mean(mean((loc_N_stst - loc_avg_n).^2)))]

equilibrium_n = [equilibrium_n, loc_avg_n]



d_rho = [d_rho, fluct_n/l^2]

eq_d_rho = [eq_d_rho, sqrt(loc_avg_n)/l^2]
  

end

 d_rho_vs_Cov_P = [d_rho_vs_Cov_P, d_rho(th_i)] 
 d_rho_eq_vs_Cov_P = [d_rho_eq_vs_Cov_P, eq_d_rho(th_i)];
 
end

cov_axis = str2double(cov);

figure(1)
h = loglog(cov_axis, d_rho_vs_Cov_P, 'o-');
legend_handles = [legend_handles, h];
hold on
loglog(cov_axis, d_rho_eq_vs_Cov_P, 'r--')

end

h = loglog(cov_axis,d_rho_eq_vs_Cov_P, 'r--');
legend_handles = [legend_handles, h];
legend_str{P_i + 1} = ['\langle\Delta{\rho}\rangle_{eq}'];

  

title('\langle\Delta{\rho}\rangle vs. coverage fraction, apolar movement')
xlabel('\kappa (coverage fraction)');
ylabel('\langle\Delta{\rho}\rangle');

ax = gca;
set(ax, 'FontSize', 16)

legend(legend_handles, legend_str,'Location', 'NorthWest')

drawnow


end




















