//
//  main.cpp
//  PA_18_01_16
//
//  Created by 11678505 on 18/01/2016.
//  Copyright (c) 2016 11678505. All rights reserved.

// this code was written by Cameron Zachreson (11678505 was my UTS student ID number)
/* This script simulates surface-motile bacteria that secrete EPS trails and deform the local topography
the code outputs a matrix of position, orientation, and length data for N bacteria in the following format:
X1(t = 1), Y1(t = 1); theta1(t = 1); length1(t = 1); X1(t = 2); Y1(t = 2); theta1(t = 2); length1(t = 2);... length1(tf)
X2(t = 1), Y2(t = 1); theta2(t = 1); length2(t = 1); X2(t = 2); Y2(t = 2); theta2(t = 2); length2(t = 2);... length2(tf)
.
.
.
XN(t = 1), YN(t = 1); thetaN(t = 1); lengthN(t = 1); XN(t = 2); YN(t = 2); thetaN(t = 2); lengthN(t = 2);... lengthN(tf)
this matrix can be read by the XYTL class for analysis in the .m files provided.
the script also creates text files containing the values of each pixel in the stigmergy grid for both the EPS and topography
the frequency of data recording can be set by changing d_rec
*/


//#include <mkl.h>
#include <cstdlib>
#include <algorithm>
#include <random>

#include <iostream>
#include <sys/stat.h>
#include <sstream>

#include <array>

#include "/home/cjzachre/Documents/headers/general_functions.h"
#include "/home/cjzachre/Documents/headers/bins_class.h"
#include "/home/cjzachre/Documents/headers/stigmergy_grid_class.h"
#include "/home/cjzachre/Documents/headers/particle_class.h"
#include "/home/cjzachre/Documents/headers/image_matrix_class.h"
#include "/home/cjzachre/Documents/headers/output_matrix_class.h"

int main(int argc, const char * argv[]) {
    
    string output_directory = "/home/cjzachre/Documents/PA_PB_followup_25_03_2016/output/";
    mkdir(output_directory.c_str(), 0777);
    
    double seed = 1;
    
    default_random_engine generator;
    generator.seed(seed); // for gaussian distribution (reversal periods)
    
    srand(seed);//for uniform distributions
    
    double s_SG = 0.25;
    
    // top level parameter is mean reversal period, for each value, iterate through gamma_s values
    double mean_reversal_period = 1000;
    double stdev_reversal_period = mean_reversal_period / 5;
    
    ostringstream label_lvl_1;
    label_lvl_1 << "T_rev_" << mean_reversal_period;
    string identifier_1 = label_lvl_1.str(); // name of the created folder
    
    string output_directory_lvl_1 = output_directory + identifier_1 + "/";
    mkdir(output_directory_lvl_1.c_str(), 0777);
    
    
    double gamma_s_coeffs[] = {1.6}; // modulates the topographical force
    
    int num_params = 1;
    
    for (int param = 0; param < num_params; param++)
    {
    
    double gamma_s_coeff = gamma_s_coeffs[param];
    
    ostringstream label_lvl_2;
    label_lvl_2 << "test_g_s_" << gamma_s_coeff;
    string identifier_2 = label_lvl_2.str(); // name of the created folder
    string output_directory_lvl_2 = output_directory_lvl_1 + identifier_2 + "/";
    mkdir(output_directory_lvl_2.c_str(), 0777);
        
    double gamma_s = gamma_s_coeff / s_SG;
    double gamma_l = 0.0 / s_SG; // gamma_l is another environmental foce like gamma_s, 
	    			 // but k_l and beta_l will be faster, (l stands for 'liquid')
	    			 // gamma_l > 0 implements cell aggregation, this was not included in the thesis.
        
    double k_p = 0.1 * (s_SG * s_SG); //deposition rate per unit area
    double k_l = 10 * (s_SG * s_SG); //doesn't do anything unless gamma_l > 0
    double k_s = 0.05 * (s_SG * s_SG) * (1 / gamma_s_coeff);
        
    double beta_p = 0.0005; //exponential degradation constant
    double beta_s = 0.00025 * (1 / gamma_s_coeff); // degradation constant for topographical trace
    double beta_l = 0.1;
    
    double c_max = 1.0 * (s_SG * s_SG);
    double r0 = 1;
    
    double Rep_const = 1;
    
    double aa = 1 / ( pow((13.0/7.0), (7.0/6.0)) * r0 );
    double bb = 1 / ( pow((13.0/7.0), (13.0/6.0)) * r0 );
    
    
    double Eo = Rep_const / (12*(aa - bb));
    
    double F_max = 10;
    double Mu = 1;
    
    double F_p = 1.5;
    
    double l_min = 3.0;
    double w = 1;
    double l_max = 2 * l_min + w;
    double s_CL = l_max + w;
    
    double phi = 0.5 * M_PI;
    double l_pil = l_min;
    double retraction_period = 10;
    
    normal_distribution<double> T_rev_dist(mean_reversal_period, stdev_reversal_period);
    
    double T_rev = T_rev_dist(generator);// a specific reversal clock which is selected from the distribution at random
    double T_rev_max = mean_reversal_period + 4 * stdev_reversal_period;
    double T_rev_min = max(0.0, mean_reversal_period - 4 * stdev_reversal_period);
    
    
    double P_min = 0.1;
    double P_max = 0.3;
    double P_attach_particle = 0.25;
    double L = s_CL * 20;
    double Ls_x = L;
    double Ls_y = L;
    int N = 1000;
    
    double dt_max = 0.1;
    double dt_min = 0.005;
    double dt = dt_max; // this will change to keep the maximum particle movement below dX_max or dTheta_max
    double t = 0.0;
    double tf = 50000;
    double big_num = 100000000;
    
    double dX_max = 0.1;
    double V_max = 0;
    
    double d_rec_XYTL = 20;
    double d_rec_image = 1000;
    double d_check = 1;
        
    double t_last_rec_XYTL = 0;
    double t_last_rec_image = 0;
    double t_last_check = 0;
        
    bool check;
    bool rec_XYTL;
    bool rec_image;
    int inc = -1;
   
    //save parameter list
    string params = output_directory_lvl_2 + "parameters" + identifier_1 + "_" + identifier_2 + ".txt";
    // writes a text file, saved to input directory, containing the specified string
    ofstream params_txt;
    params_txt.open(params.c_str());
    
            params_txt << "parameter" << "\t" << "value" << "\n"
        << "number of particles" << "\t" << N << "\n"
        << "mean reversal period" << "\t" << mean_reversal_period << "\n"
        << "standard deviation T_rev" << "\t" << stdev_reversal_period << "\n"
        << "gamma substratum coeff" << "\t" << gamma_s_coeff << "\n"
        << "gamma substratum" << "\t" << gamma_s << "\n"
        << "gamma liquid" << "\t" << gamma_l << "\n"
        << "k attachment bias" << "\t" << k_p << "\n"
        << "k substratum" << "\t" << k_s << "\n"
        << "k liquid" << "\t" << k_l << "\n"
        << "beta attachment bias" << "\t" << beta_p << "\n"
        << "beta substratum" << "\t" << beta_s << "\n"
        << "beta liquid" << "\t" << beta_l << "\n"
        << "retraction period" << "\t" << retraction_period << "\n"
        << "minimum attachment probability" << "\t" << P_min << "\n"
        << "maximum attachment probability" << "\t" << P_max << "\n"
        << "pilus-particle attachment probability" << "\t" << P_attach_particle << "\n"
        << "rng_seed" << "\t" << seed << "\n"
        << "minimum length" << "\t" << l_min << "\n"
        << "width" << "\t" << w << "\n"
        << "maximum length" << "\t" << l_max << "\n"
        << "side length sorting" << "\t" << s_CL << "\n"
        << "side length environment" << "\t" << s_SG << "\n"
        << "maximum count" << "\t" << c_max << "\n"
        << "repulsion radius" << "\t" << r0 << "\n"
        << "repulsion constant" << "\t" << Rep_const << "\n"
        << "aa(LJ)" << "\t" << aa << "\n"
        << "bb(LJ)" << "\t" << bb << "\n"
        << "Eo(LJ)" << "\t" << Eo << "\n"
        << "max repulsion force" << "\t" << F_max << "\n"
        << "friction coefficient" << "\t" << Mu << "\n"
        << "propulsion force" << "\t" << F_p << "\n"
        << "phi (pili angle)" << "\t" << phi << "\n"
        << "pilus length" << "\t" << l_pil << "\n"
        << "system side length" << "\t" << L << "\n"
        << "x dim of source area" << "\t" << Ls_x << "\n"
        << "y dim of source area" << "\t" << Ls_y << "\n"
        << "max possible time step" << "\t" << dt_max << "\n"
        << "min possible time step" << "\t" << dt_min << "\n"
        << "t_final" << "\t" << tf << "\n"
        << "dX max" << "\t" << dX_max << "\n"
        << "config recording period" << "\t" << d_rec_XYTL << "\n"
        << "image recording period" << "\t" << d_rec_image << "\n"
        << "list check period" << "\t" << d_check << "\n";
    
    params_txt.close();
    
    string XYTL_filename = output_directory_lvl_2 + "XYTL_" + identifier_1 + "_" + identifier_2 + ".txt";
    
    //create output file for images
    string dirname_s = output_directory_lvl_2 + "images_" + identifier_1 + "_" + identifier_2 + "_S/";
    mkdir(dirname_s.c_str(), 0777);
    
    string dirname_l = output_directory_lvl_2 + "images_" + identifier_1 + "_" + identifier_2 + "_L/";
    mkdir(dirname_l.c_str(), 0777);
    
    string dirname_p = output_directory_lvl_2 + "images_" + identifier_1 + "_" + identifier_2 + "_P/";
    mkdir(dirname_p.c_str(), 0777);
    
    
    //initialize output vectors for X, Y, theta
    
    output_matrix X_Y_T_L(tf, d_rec_XYTL, N, 4);
    
    // initialize bins
    BINS bins(L, L, s_CL);
    
    //initialize stigmergy grid
    Stigmergy_Grid stig_grid(L, L, s_SG, w, P_min, P_max, gamma_s, gamma_l);
    
    // initialize N particles
    vector <Particle> particles;
    
    for (int i = 0; i < N; i++)
    {
        double x = ( rng() - 0.5 ) * (Ls_x);
        double y = ( rng() - 0.5 ) * (Ls_y);
        double theta = rng() * 2 * M_PI;
	
        double l = l_min + rng() * (l_max - l_min);
        
        T_rev = rng() * mean_reversal_period;
        while (T_rev < T_rev_min || T_rev > T_rev_max)
        {T_rev = T_rev_dist(generator);}
        
        Particle p(x, y, l, phi, theta, retraction_period, T_rev, l_pil, F_p, i,
                   N, bins.XI_CL_max, bins.YI_CL_max, L, L);
        
        p.bin(s_CL, L, L);
        
        p.stig_index(s_SG, L, L);
        
        stig_grid.index_particle_to_grid(p.xi_sg, p.yi_sg, p.theta, p.l, w, s_SG, L, L);
        
        stig_grid.add_counts(k_p, k_s, k_l, dt, c_max);
        
        bins.put_in_bin(p.YI, p.XI, p.ID);
        
        particles.push_back(p);
    }
    
    for (int a = 0; a < mean_reversal_period * 10; a++)
    {
        for (int i = 0; i < N; i++)
        {
            
            T_rev = T_rev_dist(generator);
            while (T_rev < T_rev_min || T_rev > T_rev_max)
            {T_rev = T_rev_dist(generator);}
            
            particles[i].init_reversal_dist(T_rev);
            
        }
    }
    
        
    for (int step = 0; step < big_num; step++)
    {
        
        check = 0;
        rec_XYTL = 0;
	rec_image = 0;
        V_max = 0;
        
        t += dt;
        
        if (t >= tf){break;}
        
        if (t - t_last_check > d_check){t_last_check = round(t); check = 1;}
        
        if (t - t_last_rec_XYTL > d_rec_XYTL)
        {
	    cout << setprecision(9)<< t << "\t" << dt << "\t" << step << endl;
            t_last_rec_XYTL = round(t); rec_XYTL = 1; inc += 1;
            cout << setprecision(9)<< t << "\t" << step << endl;
        }
            
        if (t - t_last_rec_image > d_rec_image)
        {
            t_last_rec_image = round(t); rec_image = 1;
        }
        
        
        
        // reset order parameter for this timestep
        double cos_t_tot = 0;
        double sin_t_tot = 0;
        
        stig_grid.decay(beta_p, beta_s, beta_l, dt);
        
        for (int i = 0; i < N; i++)
        {
	
	   
            particles[i].move_and_reset(dt, N, L, L);
         
	   
	   
            particles[i].bin(s_CL, L, L);
	    
	    
	    
            bins.put_in_bin(particles[i].YI, particles[i].XI, particles[i].ID);
	    
            //update stigmergy grid for the current particle
            particles[i].stig_index(s_SG, L, L);
	   
            stig_grid.index_particle_to_grid(particles[i].xi_sg, particles[i].yi_sg,
                                             particles[i].theta, particles[i].l, w, s_SG, L, L);
          
            particles[i].F_x += round(stig_grid.F_env_x * 100000) / 100000;
            particles[i].F_y += round(stig_grid.F_env_y * 100000) / 100000 ;
            particles[i].torque += round(stig_grid.torque_env * 100000) / 100000;
   
            stig_grid.clear_for_next_particle();
            
            cos_t_tot += cos(particles[i].theta);
            sin_t_tot += sin(particles[i].theta);
        }
        
        
        stig_grid.add_counts(k_p, k_s, k_l, dt, c_max);
        
        stig_grid.update_stig_grid();

        for (int i = 0; i < N; i++)
        {
            // store data for this timestep
            //cout << i << endl;
            if (rec_XYTL)
                
            {
                X_Y_T_L.insert_data(inc, i, 0, particles[i].X);
                X_Y_T_L.insert_data(inc, i, 1, particles[i].Y);
                X_Y_T_L.insert_data(inc, i, 2, particles[i].theta);
                X_Y_T_L.insert_data(inc, i, 3, particles[i].l);
            }
            
            particles[i].Repulsion(particles, bins, i, r0, Eo, F_max, L, L, check);
            
            
            T_rev = T_rev_dist(generator);
            while (T_rev < T_rev_min || T_rev > T_rev_max)
            {
                T_rev = T_rev_dist(generator);
            }
            
            particles[i].Twitch(dt, L, L, s_SG, stig_grid, s_CL, bins, particles, P_attach_particle, T_rev);
            
        }
        
        
        bins.clear_bins();
        
        for (int i = 0; i < N; i++)
        {
            particles[i].V_tot(Mu);
            V_max = max(V_max, particles[i].V_pole_max);
        }
        
        dt = min(dX_max/V_max, dt_max);
        dt = max(dt, dt_min);
        
        if (rec_image)
        {
            
            string time = ZeroPadNumber((int)round(t), 5);
            
            cout << time << endl;
            
            string filename_l = dirname_l + time + ".txt";
            string filename_p = dirname_p + time + ".txt";
            string filename_s = dirname_s + time + ".txt";
            
            image_matrix stig_map_L(filename_l, stig_grid.Stig_Grid_L_now, stig_grid.XI_SG_max, stig_grid.YI_SG_max);
            image_matrix stig_map_P(filename_p, stig_grid.Stig_Grid_P_now, stig_grid.XI_SG_max, stig_grid.YI_SG_max);
            image_matrix stig_map_S(filename_s, stig_grid.Stig_Grid_S_now, stig_grid.XI_SG_max, stig_grid.YI_SG_max);
        }
	
        
    }
        
    //write data to text file
    
    X_Y_T_L.write_to_file(XYTL_filename);
    
    cout << L/2 << endl;
        
    }
    
    return 0;
}
