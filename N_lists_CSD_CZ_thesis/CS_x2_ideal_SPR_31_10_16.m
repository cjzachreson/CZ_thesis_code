% this code was used to find cluster size distributions
% for variants of the twitching rod model. The data is
% presented in Fig. 4.9 of the thesis.

close all
clear all


P = {'1'}
cov = {'0.3'};
intervals = {'20000','30000','40000','50000','60000',...
             '70000','80000', '90000', '100000'};

modes = {'SPR'}%'SPR', 'Fp_1p5', 'len_dist', 'phi', 'Ppp_0p25', 'P_a_0p1'}%{'apolar}

for m_i = 1:length(modes)
    mode = modes{m_i};
    

for P_i = 1:length(P)
    
    legend_str = cell(1, numel(cov));
    legend_handles = [];
    
for cov_i = 1:length(cov)

    
    if strcmp(mode, 'P_a_0p1') == 1
        P{P_i} = '0.1';
    else
        P{P_i} = '1';
    end
    


['P = ' P{P_i}]
['cov = ' cov{cov_i}]



root = ['/Users/11678505/Desktop/live/simulations/TR_PB_no_stig/'...
    'TR_PB_SPR_model_tests_output/AR_6_kappa_0p3_' modes{m_i} '/AR_6_kappa_0p3_'...
    modes{m_i} '_L_240_P_' P{P_i}];
           
        
    output_path = [root '/x2_bin_CS_gl_cov_0p3_P_' P{P_i} '_stst/'];
    if ~isdir(output_path)
    mkdir(output_path);
    end
    
    CS_Mat = [];
    
for int_i = 1:length(intervals)
 
    intervals{int_i}
    
    lists_path = [root ...
    '/N_list_cut_1.5_trep_' intervals{int_i} '_cov_' cov{cov_i}...
    '_P_' P{P_i} '/'] ;  
    
    
    N_lists = dir(lists_path);
    d_rec = 20;
    
    
    
    
file = zeros(1, length(N_lists));

for i = 1:length(N_lists)
    
    try
        check_1 = N_lists(i).name(end - 3: end) == '.txt';
        check_2 = N_lists(i).bytes > 0;
        check = check_1 & check_2;
        if check
            file(i) = 1;
        end
    end
    
end

N_lists = N_lists(file == 1);
    
    
for entry = 1:length(N_lists) 
        
      t = entry * d_rec;
       
      N_list = dlmread([lists_path '/' N_lists(entry).name]);
      N_list = N_list';
      N = length(N_list);
   
        
%below is the cluster analysis code

Call = zeros(1, N);
ClusterSize = zeros(1, N);

    for i = 1 : N
    
        if Call(i) == 0

        j = nonzeros(N_list(:, i));

            if full(j)

                jj = sort(nonzeros(N_list(:, j)));

                jj = jj(diff([0 ; jj])~=0);

                if full(jj)

                while ~isequal(j, jj)

                     j = [j ; jj];

                     j = sort(nonzeros(j));
                     j = j(diff([0 ; j])~=0);


                     jj = sort(nonzeros(N_list(:, j)));
                     jj = jj(diff([0 ; jj])~=0);


                end

                end

            ClusterSize(jj) = numel(jj);

            Call(jj) = 1;


            else 
                ClusterSize(i) = 1;
            end


        end
    
    end
    
    lastedge = 1;
    edges = 1;
    inc = 1;
    while edges(end) < N
        
       if inc > 1 
       edges = [edges, inc^2];
       end
       
       inc = inc+1;
       
       
    end
    
    
    
    CS_mat_t = hist(ClusterSize, edges); 
    
    CS_Mat = [CS_Mat ; CS_mat_t./[1, diff(edges)]/N];
    
   
    
end
    
end

%   imagesc(CS_Mat)

  
     CS_name = [output_path 'netCS.txt'];
     dlmwrite(CS_name, CS_Mat)
     
     bin_vals_name = [output_path 'x_bin_edges.txt'];
     dlmwrite(bin_vals_name, edges)
     
     bin_width_name = [output_path 'x_bin_width.txt'];
     dlmwrite(bin_width_name, [1, diff(edges)])
    
     
     
    t_integrated_CS = sum(CS_Mat)./500;
    
    
    
    legend_str{cov_i} = ['P = ' P{P_i} ', ' '\rho = ' cov{cov_i}];
    
    a = semilogy(edges, t_integrated_CS, 'o-');
    legend_handles = [legend_handles, a];
    hold on
%     limits = [0, edges(end), 0, 1];
%     axis(limits)
    title(['CSD, ' mode])
    xlabel('cluster size');
    ylabel('P(CS)');
    legend(legend_handles, legend_str{1:cov_i}) 
    ax = gca;
    set(ax, 'FontSize', 16)
    drawnow 
    
     saveas(a, [output_path 'integrated_CS_logy.pdf'], 'pdf')
%     
     integrated_CS_name = [output_path 'integrated_CS.txt'];
     dlmwrite(integrated_CS_name, t_integrated_CS)
    
    
  

hold off

end
end
end
    
       
    
        
    
    
        

