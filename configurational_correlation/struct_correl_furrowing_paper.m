%structural correlation takes a coarse grid of all the occupied space in a
%frame and compares it to the subsequent frames. This gives a decay of
%structural information with time that can be used as a direct measure of
%trail following behavior. Must be compared with a similar metric that
%tracks whether or not specific particles have changed their coarse
%locations for unambiguous quantification of path following. IE -
%structural decay time much longer than particle residence time. 

close all
clear all

gl = {'0'} 
state = { 'WT' };
Trev = {'1000'};
%gs = {'0.001','0.25','0.5','0.6','0.7','0.8','1','1.1','1.2','1.3','1.4','1.5','1.6'};
%gs = {'0.9','1','1.1','1.3'};
gs = {'2'}
listindex = 0;
for gs_i = 1:length(gs)
for gl_i = 1:length(gl)
for s_i = 1:length(state)
for tr_i = 1:length(Trev)

listindex = listindex + 1;
    
% root = ['/Users/11678505/Desktop/PA_sim_data/PB_followup/tf_100k/'];
root = ['/Users/11678505/Desktop/PA_sim_data/PB_followup/'];
XYTL_filename = [root ...
    'XYTL_T_rev_' Trev{tr_i} '_test_g_s_' gs{gs_i} '.txt'];

output_path = [root 'decay_constants/']
if ~isdir(output_path)
    mkdir(output_path)
end

L = 160

XYTL = XYTL_data(XYTL_filename, L/2);

dims = size(XYTL.l);

numframes = dims(2) - 1;

N = dims(1);

grid_res = L / 40;

XI_max = L / grid_res;
YI_max = L / grid_res;

structure_grid = zeros(YI_max, XI_max, numframes);
exchange_grid = zeros(N, numframes);

XI = floor(XYTL.X / grid_res) + (L / grid_res / 2) + 1;
YI = abs(ceil(XYTL.Y / grid_res) - (L / grid_res / 2)) + 1;

XI(XYTL.X == L/2) = XI_max;
YI(XYTL.Y == - L/2) = YI_max;

exp_fit = fittype('exp_decay(x, b)');
%dbl_diff_exp_fit = fittype('dbl_diff_exp_decay(x, b, c)');

%fill structure grid
for t = 1:numframes
    
    for i = 1:N
        structure_grid(YI(i, t), XI(i, t), t) = 1;
    end
    
end

%fill exchange grid
for t = 2:numframes
    for i = 1:N
        if XI(i, t) == XI(i, t - 1) && YI(i, t) == YI(i, t-1) 
            exchange_grid(i, t) = 1;
        end
    end
end
                   

%calculate structural correlation values and positional correlation values
structural_correlation = zeros(numframes, numframes);
positional_correlation = zeros(numframes, numframes);
decay_constants_struct = [];
decay_constants_pos = [];
decay_constants_ratio = [];

for t = 1900:numframes-100
    t
    current_structure = structure_grid(:, :, t);
    current_positions = exchange_grid(:, t);
    
    correl_structure_back = current_structure;
    correl_structure_forward = current_structure;
    
    correl_positions_back = current_positions;
    correl_positions_forward = current_positions;
    
    n_pos_ti = sum(current_positions);
    n_struct_ti = sum(sum(current_structure));
    
    inc = 0;
    
    for t_then = 1 : t-1
        inc = inc + 1;
           
            correl_structure_back = correl_structure_back .* structure_grid(:, :, t - t_then);
            structural_correlation(t, (t - inc)) = sum(sum(correl_structure_back)) / n_struct_ti;
            
            correl_positions_back = correl_positions_back .* exchange_grid(:, t - t_then);
            positional_correlation(t, (t - inc)) = sum(correl_positions_back) / n_pos_ti;
            
    end
    
    for t_next = t : numframes
        inc = inc + 1;
            
            correl_structure_forward = correl_structure_forward .* structure_grid(:, :, t_next);
            structural_correlation(t, inc) = sum(sum(correl_structure_forward)) / n_struct_ti;
            
            correl_positions_forward = correl_positions_forward .* exchange_grid(:, t_next);
            positional_correlation(t, inc) = sum(correl_positions_forward) / n_pos_ti;
    end
    f1 = structural_correlation(t, max([1, t-100]):min([t+100, numframes]));  
    f2 = positional_correlation(t, max([1, t-100]):min([t+100, numframes]));
    
    f1_forward = structural_correlation(t, t : min([t+100, numframes]));
    x1 = 0:(numel(f1_forward))-1;
    
    f2_forward = positional_correlation(t, t : min([t+100, numframes]));
    x2 = 0:(numel(f1_forward))-1;
    
    f12_forward = f1_forward - f2_forward;
    x12 = 0:numel(f12_forward) - 1;
    
     fit_f1 = fit(x1', f1_forward', exp_fit, 'StartPoint', -0.25 );
     b1 = coeffvalues(fit_f1);
     decay_constants_struct = [decay_constants_struct, b1];
%    
     fit_f2 = fit(x2', f2_forward', exp_fit, 'StartPoint', -0.25 );
     b2 = coeffvalues(fit_f2);
     decay_constants_pos = [decay_constants_pos, b2];

     decay_constants_ratio = [decay_constants_ratio, b2/b1];
    
     fit_12 = exp(b1 .* x12) - exp(b2.*x12);
     
%     figure(50)
%     plot(x12, f12_forward, '--k')
%     hold on
%     plot(x12, fit_12, 'r')
%     hold off
    

end
%figure(listindex)
% edge_min = mean(decay_constants_struct) - 5 * std(decay_constants_struct);
% edge_max = mean(decay_constants_struct) + 5 * std(decay_constants_struct);
% n_bins = 100;
% edges = [edge_min:((edge_max - edge_min)/n_bins):edge_max];
% histogram(decay_constants_struct, edges)

% edge_min = mean(decay_constants_ratio) - 5 * std(decay_constants_ratio);
% edge_max = mean(decay_constants_ratio) + 5 * std(decay_constants_ratio);
% n_bins = 100;
% edges = [edge_min:((edge_max - edge_min)/n_bins):edge_max];
% histogram(decay_constants_ratio, edges)

drawnow
dlmwrite([output_path 'decay_constants_struct_gs_' gs{gs_i} '.txt'], decay_constants_struct)
dlmwrite([output_path 'decay_constants_pos_gs_' gs{gs_i} '.txt'], decay_constants_pos)
%dlmwrite([output_path 'decay_constants_ratio_gs_' gs{gs_i} '.txt'], decay_constants_ratio)
end
end
end
end












