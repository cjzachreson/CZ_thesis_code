% path following - stigmergy. Stigmergy is an individual modifying the
% local environment to influence the behavior of the same or a different
% individual arriving at the SAME area at a LATER point in time. This code
% will analyze the movement vectors of particles that enter into local
% regions and compare the entry vector to the vector of the last particle
% that entered the area. The analysis is only done upon particle ENTRY into
% a local area, thus ignoring the effects of dwell. The dot products of
% these unit vectors should approach a time average of 0 if no
% trail-following or drift velocity exisists. 

close all
clear all

Pmin = {'0.1'}%, '1'};
Pmax = {'1'};
Ppp = {'0'}

rep = {'repulsion'}%{'no_repulsion', 'repulsion'}

tret = {'10'}%, '10'}
Trev = {'1000'};
kpfac = {'0.1', '0.05', '0.01', '0.005', '0.001'};
b_p = {'0.01', '0.005', '0.001', '0.0005', '0.0001'};
Ns = {'250'}%, '250', '500'};

tfinal = '50k'

index = 0;

output_index = 0

delay_periods = [10, 20, 40, 80, 160, 320, 640];

for dt_i = 1:numel(delay_periods)

    output_index = output_index + 1;
    
for rep_i = 1:numel(rep)
for pmin_i = 1:numel(Pmin)
for pmax_i = 1:numel(Pmax)
for pp_i = 1:numel(Ppp)
    
Pstate =  ['Pmin_' Pmin{pmin_i} '_Pmax_' Pmax{pmax_i} '_Ppp_' Ppp{pp_i} ] ;

for tr_i = 1:numel(tret)
for N_i = 1:length(Ns)
    
    
   pow_PF = NaN(numel(kpfac), numel(b_p));

   y0_PF = NaN(numel(kpfac), numel(b_p));

   Ao_PF = NaN(numel(kpfac), numel(b_p));
    
    
    row_dc = 0;
    
    x = [];
    y = [];
    v1 = [];
    v2 = [];
    
for k_i = 1:length(kpfac)
    
    row_dc = row_dc + 1;
    
    col_dc = 0;
    
for bp_i = 1:length(b_p)
    
    col_dc = col_dc + 1;
    
root = ['/Users/11678505/Desktop/live/simulations/EPS_only/output/'...
    'tf_' tfinal '/' rep{rep_i}...
    '/Trev_' Trev{1} '_tret_' tret{tr_i} '_N_' Ns{N_i} '/' ...
    Pstate '_bp_' b_p{bp_i} '_kpfac_' kpfac{k_i} '/'];


if isdir(root)

 index = index + 1

XYTL_filename = [root ...
    'XYTL_Trev_1000_tret_' tret{tr_i} '_N_' Ns{N_i} '_' Pstate ...
    '_bp_' b_p{bp_i} '_kpfac_' kpfac{k_i} '.txt'];

output_path = [root 'path_following/'];
if ~isdir(output_path)
    mkdir(output_path)
end

d_rec = 20;

L = 112;

XYTL = XYTL_data(XYTL_filename, L/2);

dims = size(XYTL.l);

numframes = dims(2) - 1;

t_space = 10;

interval = round(numframes/1.5) : t_space : numframes;

%delay_periods = 700%100;%300:100:(numel(interval) * t_space);%[1, 10:20:80, 100:100:(numel(interval) * t_space)];

N = dims(1);

IDs = 1:N;

grid_res = L / 28;

XI_max = L / grid_res;
YI_max = L / grid_res;

% orientation_grid = zeros(YI_max, XI_max, numframes);
orientation_grid = NaN(YI_max, XI_max, numframes);
% records orientations of the last particle to pass through each space

exchange_grid = zeros(N, numframes);

XI = floor(XYTL.X / grid_res) + (L / grid_res / 2) + 1;
YI = abs(ceil(XYTL.Y / grid_res) - (L / grid_res / 2)) + 1;

XI(XYTL.X == L/2) = XI_max;
YI(XYTL.Y == - L/2) = YI_max;

% fill entry grid, if a particle enters a new region it gets a 1, if it
% stays put it gets a 0

exchange_grid(:, 1) = 0; % all particles enter at initialization

for t = 2:numframes
    for i = 1:N
        if XI(i, t) ~= XI(i, t - 1) || YI(i, t) ~= YI(i, t-1) 
            exchange_grid(i, t) = 1;
        end
    end
end

% fill orientation grid
for t = 1:numframes
    for i = IDs(exchange_grid(:, t) == 1)
        orientation_grid(YI(i, t-1), XI(i, t-1), t-1) = XYTL.theta(i, t-1);
    end
end

% want to look at the orientational memory in the system, 
% but only upon particle entry, ie, ignoring dwell   
% when a particle changes its index, compare its orientation to the
% orientation of the local memory vector (the orientation of the last
% particle that entered) 

PF_vs_dt = NaN(1, numel(interval));

exchange_frac_vs_dt = NaN(1, numel(interval));

exchanges_vs_dt = NaN(1, numel(interval));

PF_index = 0;

for t = interval
    
    PF_index = PF_index + 1;
    
    delay = delay_periods(dt_i);
        
    exchanges_dt = 0;

    exchange_frac_dt = 0;

    if t + delay < numframes

        theta_grid_t_dt = orientation_grid(:, :, [t : t+delay]);

        vec_x = cos(2 * theta_grid_t_dt);
        vec_y = sin(2 * theta_grid_t_dt);
        
        exchanges_dt = sum(sum(exchange_grid(:, [t : t+delay])));

        exchange_frac_dt = ...
            exchanges_dt ./ (N * numel(t : t+delay));

        exclusion_mat = theta_grid_t_dt;
        exclusion_mat(isnan(exclusion_mat)) = 0;
        exclusion_mat = sum(sign(exclusion_mat), 3);
        
        exclusion_mat(exclusion_mat <= 1) = NaN;
        
        exclusion_mat = sign(exclusion_mat);
        
        avg_vec_x = nanmean(vec_x, 3) .* exclusion_mat;

        avg_vec_y = nanmean(vec_y, 3) .* exclusion_mat;

        avg_vec_mag = sqrt(avg_vec_x.^2 + avg_vec_y.^2);

        PF_vs_dt(PF_index) = nanmean(nanmean(avg_vec_mag));

        exchange_frac_vs_dt(PF_index) = exchange_frac_dt;

        exchanges_vs_dt(PF_index) = exchanges_dt;
        
    end
    
%      figure(300)
%      XYTL.Quiver(t-1, 1, 0)
%      drawnow
    
end

% average exchange coherence over each time interval
PF_avg_x_t = nanmean(PF_vs_dt);

t_vals = delay_periods(~isnan(PF_avg_x_t)) * d_rec;

PF_avg_x_t = PF_avg_x_t(~isnan(PF_avg_x_t));

x = [x, log(str2double(b_p{bp_i}))];
y = [y, log(str2double(kpfac{k_i}))];
v1 = [v1, PF_avg_x_t];

  
 
end
end
end


space = 0.01;
buffer = space;

kq = (floor(min(x)/space) * space - buffer) :...
    space : (ceil(max(x)/space) * space + buffer);

gq = (floor(min(y)/space) * space - buffer) :...
    space : (ceil(max(y)/space) * space + buffer);

[xq, yq] = meshgrid(kq, gq);

vq1 = griddata(x, y, v1, xq, yq, 'natural');%'natural','v4','cubic'


figure(output_index)

contourf(xq, yq, vq1)%, 'LevelStep', abs(min(min(vq))/10))
colorbar
caxis([0 1])
hold on
plot(x, y, 'o', 'color', 'r')
title(['time interval = ' num2str(delay * d_rec) 's'])
drawnow


end
end
end
end
end
end

    


end





