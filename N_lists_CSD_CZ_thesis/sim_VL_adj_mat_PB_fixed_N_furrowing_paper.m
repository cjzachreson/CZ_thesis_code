% 
% This script was used to find neighbor lists
% for cluster analysis, this specific version
% of the code was used for some of the data
% presented in Chapter 4.4 of the thesis

close all
clear all

gl = {'0'} 
state = { 'WT' }
Trev = {'1000'}
gs = {'0.9', '1', '1.1', '1.3'}
listindex = 0;

for gl_i = 1:length(gl)
for s_i = 1:length(state)
for tr_i = 1:length(Trev)
for gs_i = 1:length(gs)
listindex = listindex + 1



threshold = '1.5'

root = ['/Users/11678505/Desktop/PA_sim_data/PB_followup/tf_100k/'];

N_list_folder_name = [root ...
    'N_list_cut_' threshold '_gl_' gl{gl_i}...
    '_TR_' Trev{tr_i} '_gs_' gs{gs_i} '_' state{s_i} '/'] ;    
        
XYTL_filename = [root ...
    'XYTL_T_rev_' Trev{tr_i} '_test_g_s_' gs{gs_i} '.txt'];

if ~isdir(N_list_folder_name)
mkdir(N_list_folder_name)
end

L = 80;
XYTL = XYTL_data(XYTL_filename, L);


dims = size(XYTL.l);
tf = dims(2)
Lx = 160
Ly = 160
d_rec = 20;

%Cell lists (neighborhood grid)
s_CL = 8;%side length - set to maximum length of a cell
XI_CL_max = Lx / s_CL;
XI_CL_min = 1;
YI_CL_max = Ly / s_CL;
YI_CL_min = 1;


XI_CL = ceil( XYTL.X / s_CL ) + ceil(Lx / s_CL / 2);
YI_CL = ceil(Ly / s_CL / 2) - floor(XYTL.Y / s_CL);



for t = 1:tf
    
    if rem (t * d_rec , 1000) == 0
    t * d_rec
    end
    
    N = nnz(XYTL.l(:, t));
    IDs = 1:N;
    
    XI_CL_t = ceil( XYTL.X(IDs, t) / s_CL ) + ceil(Lx / s_CL / 2);
    
    YI_CL_t = ceil(Ly / s_CL / 2) - floor(XYTL.Y(IDs, t) / s_CL);
    
    N_list_t = zeros(N, 10);
    
    
    for i = IDs
       
        N_list_i_t = zeros(1, 10);
       
        hood = neighborhood_PB(IDs, XI_CL_t', YI_CL_t', ...
                       XI_CL_max, XI_CL_min, YI_CL_max, YI_CL_min, i);
        inc = 0;
        
        for j = hood
            
            
           [ dXP, dYP ] = ...
               periodic_bounds( XYTL.X(i, t), XYTL.X(j, t), XYTL.Y(i, t), XYTL.Y(j, t), Lx, Ly );
           
           d_ijC = sqrt(dXP^2 + dYP^2);
           
           if d_ijC < 8
                  
           dist_ji = VL_Rod_dist_PB_growth(dXP, dYP, XYTL.theta(i, t),...
                     XYTL.theta(j, t), XYTL.l(i, t), XYTL.l(j, t));
          
               if dist_ji < str2double(threshold)
              
                   inc = inc + 1  ;  
              
                   N_list_i_t(inc) = j;

               end
       
           end

        end
        
    N_list_t(i, (1:numel(N_list_i_t))) = N_list_i_t; 
    end
    
    pad = zeros(1, numel(num2str(tf)) - numel(num2str(t)));
 
    list_mat_name = [N_list_folder_name, '_t_' num2str(pad) num2str(t) '.txt'];
    dlmwrite(list_mat_name, N_list_t)
    
end

end

end

end

end
