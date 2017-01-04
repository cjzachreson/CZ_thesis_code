% numerical determination of diffusion coefficient for off-lattice random
% walk
%function twitching_rods_rand_wlk_D_vs_P

clear all
close all


numreps = 1000;
numbins = numreps/10;

S_uni = log(numbins);
dS_min = (numbins/numreps) * ( 0.37 * numbins^(-0.38) + 0.045 );
dS_tau_criteria = dS_min + dS_min * 0.05;

X_i = zeros(numreps, 1);
Y_i = zeros(numreps, 1);
theta_i = ones(numreps, 1) .* 0;

ti = 0;
dt = 0.25;
tf = 100000;
dt_dS = 10;

 D_fit = fittype('D_vs_t(x, D)');

%free parameters
P_attchs = [0.05:0.05:1];
T_revs = [250, 500, 1250, 2500, 0];

%fixed parameters
l = 5;
w = 1;
F_p = 1.5;
mu = 1;
T_ret = 10;
l_pil = 5;
phi = pi/2;
attch = 0;

%three derived quantities
D_x_mat = zeros(numel(P_attchs), numel(T_revs));
D_r_mat = zeros(numel(P_attchs), numel(T_revs));
tau_mat = zeros(size(D_x_mat));


fig_index = 1;


for T_rev_i = 1:numel(T_revs)
    
    for P_attch_i = 1:numel(P_attchs)


        T_rev = T_revs(T_rev_i);

        reversals = 1;

        if T_rev == 0
            reversals = 0;
        end

        P_attch = P_attchs(P_attch_i);


        % root = ['/Users/11678505/Desktop/PA_sim_data/twitching_rods'...
        % '/Diffusion/disp_vs_t/'];
        root = ['/home/cjzachre/Documents/twitching_rods_output/Diffusion'...
            '/D_vs_P_dt_' num2str(dt) '_n_' num2str(numreps) '/'];
        if ~isdir(root)
            mkdir(root)
        end

        filename = [root, 'T_rev_' num2str(T_rev) '/'];
        if ~isdir(filename)
            mkdir(filename)
        end


        numsteps = numel(ti:dt:tf);

        mean_disp = zeros(1, numsteps);
        mean_disp_theta = zeros(1, numsteps);


        inc_dS = 0;
        dS = zeros(1, numel(ti:dt_dS:tf));

        %mean_theta_norm = VoM;


        X = X_i;
        Y = Y_i;
        theta = theta_i;

        retraction_clock = zeros(numreps, 1);
        reversal_clock = rand(numreps, 1) .* (T_rev * 2);
        r_pil = zeros(numreps, 1);
        attch = zeros(numreps, 1);

        attch_pt_x_r = zeros(numreps, 1);
        attch_pt_y_r = zeros(numreps, 1);


        for t = ti:dt:tf
            inc = round(t/dt) + 1;

            retraction_clock = retraction_clock - dt;
            reversal_clock = reversal_clock - dt;

            for i = 1:numreps

            if r_pil <= w/2
                attch(r_pil < w/2) = 0;
                retraction_clock(r_pil < w/2) = 0;
            end 

            if retraction_clock(i) <= 0

                %%%comment to disable reversals
                if reversal_clock(i) <= 0 && reversals
                    theta(i) = theta(i) + pi * sign(rand - 0.5);
                    reversal_clock(i) = rand * (T_rev * 2);
                end
                %%%

                retraction_clock(i) = rand * T_ret;
                r_pil(i) = rand * l_pil;
                phi_pil = (rand - 0.5) * phi;
                attch(i) = (rand > (1 - P_attch));

                attch_pt_x_o = (l + w) / 2 + r_pil(i) * cos(phi_pil);
                attch_pt_y_o = r_pil(i) * sin(phi_pil);

                attch_pt_x_r(i) = attch_pt_x_o * cos(theta(i)) + ...
                    attch_pt_y_o * -sin(theta(i)) + X(i); 

                attch_pt_y_r(i) = attch_pt_x_o * sin(theta(i)) + ...
                    attch_pt_y_o * cos(theta(i)) + Y(i);
            end
            end

            pole_x = (l + w) / 2 .* cos(theta) + X;
            pole_y = (l + w) / 2 .* sin(theta) + Y;

            r_pil = sqrt((pole_x - attch_pt_x_r).^2 + (pole_y - attch_pt_y_r).^2);

            twitch_vec_x = (attch_pt_x_r - pole_x) ./ r_pil;
            twitch_vec_y = (attch_pt_y_r - pole_y) ./ r_pil;

            twitch_vec_perp = cos(theta) .* twitch_vec_y - sin(theta) .* twitch_vec_x;

            F_trans_x = F_p .* twitch_vec_x;
            F_trans_y = F_p .* twitch_vec_y;
            torque = F_p .* twitch_vec_perp .* l/2;


            dX_dt = F_trans_x ./ (mu * l) .* attch;
            dY_dt = F_trans_y ./ (mu * l) .* attch;
            d_theta_dt = 12 .* torque / (mu * l^3) .* attch;

            theta = theta + d_theta_dt * dt;

            X = X + dX_dt * dt;
            Y = Y + dY_dt * dt;


            mean_disp(1, inc) = mean(sqrt((X_i - X).^2 + (Y_i - Y).^2));

            mean_disp_theta(1, inc) = mean(abs(theta_i - theta));



            if rem(t, dt_dS) == 0

                inc_dS = inc_dS + 1;

                theta_norm = atan2(sin(theta), cos(theta));


                S_theta = S_dist(theta_norm, -pi, pi, numbins);

                dS(inc_dS) =  1 - S_theta / S_uni;

            end


        end



         %find angular relaxation time
        tau = (min(find(dS <= dS_tau_criteria) ) ) * dt_dS;
        tau_mat(P_attch_i, T_rev_i) = tau;

        if T_rev ~= 0;

            t_cut = min(T_rev * 2, tau);

            else t_cut = tau;

        end


        %rotational diffusion coefficient
        f_D_r = fit([ti:dt:tf]', mean_disp_theta', D_fit,...
            'StartPoint', 0.1, 'Lower', 0.00001, 'Upper', 1000000);

        D_r = coeffvalues(f_D_r) / 2 %because <dtheta^2> = sqrt(2 * D_r * t)

        D_r_mat(P_attch_i, T_rev_i) = D_r; 


        %translational diffusion coefficient
        weights = [zeros(1, round(t_cut/dt)), ...
                   ones(1, (round(tf/dt) - round(t_cut/dt) + 1))];      

        f_D_x = fit([ti:dt:tf]', mean_disp', D_fit, ...
            'StartPoint', 1, 'Weights', weights, 'Lower', 0.00001, 'Upper', 1000000 );
        D_x = coeffvalues(f_D_x);

        D_x_mat(P_attch_i, T_rev_i) = D_x;

    end

    figure(4)
    plot(P_attchs, D_x_mat(:, T_rev_i), '*')
    saveas(gcf, [filename, 'D_x_vs_P.png'], 'png')

    figure(5)
    plot(P_attchs, D_r_mat(:, T_rev_i), '*')
    saveas(gcf, [filename, 'D_r_vs_P.png'], 'png')

    figure(6)
    plot(P_attchs, tau_mat(:, T_rev_i), '*')
    saveas(gcf, [filename, 'tau_vs_P.png'], 'png')

    %tau_name = [filename, 'tau_vs_P.txt'];
    %dlmwrite(tau_name, [tau_mat(:, T_rev_i), P_attchs'])

    %D_r_name = [filename, 'D_r_vs_P.txt'];
    %dlmwrite(D_r_name, [D_r_mat(:, T_rev_i), P_attchs'])

    %D_x_name = [filename, 'D_x_vs_P.txt'];
    %dlmwrite(D_x_name, [D_x_mat(:, T_rev_i), P_attchs'])


end

%end
  
