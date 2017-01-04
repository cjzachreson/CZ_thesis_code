
//  Created by 11678505 on 21/09/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
//

#include <cstdlib>
#include <algorithm>
#include <random>

#include <iostream>
#include <sys/stat.h>
#include <sstream>

#include <array>



#include "general_functions.h"
#include "bins_class.h"
#include "particle_class_no_stig.h"
#include "output_matrix_class.h"


int main(int argc, const char * argv[]) {
    
    if (argc != 3){
	 cout << "Usage:" << argv[0] << " P_attach L/s_CL" << endl;
	 cout << "args:" << endl;
	 for (int i = 0; i < argc; i++) 
	 { 
	  
	  cout << "i = " << i << ", arg = " << argv[i] << endl;
	  
	 }
	return 1;
	}


    double P_attach = stod(argv[1]); // this is a free parameter, the pilus-surface attachment probability
    double L_fac = stod(argv[2]); // the fraction of area covered by particles



    double l_min = 3.0;
    double w = 1;
    double l_max = 2 * l_min + w;
    double l_avg = l_min + (l_max - l_min)/2;// set initialization to use l_avg

    double kappa = l_avg + w * 2;

    double s_CL = l_max + w;




    double L = L_fac * s_CL;

    double seed = 1;
    
    string state = "AR_6_kappa_0p3_apolar";
    
    double mean_reversal_period = 1000;
    double stdev_reversal_period = 200;
    
    
    string output_directory = "/home/cjzachre/Documents/TR_PB_SPR_L_var_output/" + state + "/";
    mkdir(output_directory.c_str(), 0777);
    
    
    default_random_engine generator;
    generator.seed(seed); // for gaussian distribution (reversal period)
    srand(seed);//for uniform distributions
    
    
    
    ostringstream label_lvl_0;
    ostringstream label_lvl_1;
    
    label_lvl_0 << "_L_" << L;
    label_lvl_1 << "_P_" << P_attach;
    
    string id_0 = label_lvl_0.str();
    string id_1 = label_lvl_1.str();
    
    string output_directory_lvl_2 = output_directory + state + id_0 + id_1 + "/";
    mkdir(output_directory_lvl_2.c_str(), 0777);
    
    double r0 = 1;
    
    double Rep_const = 1;
    
    double aa = 1 / ( pow((13.0/7.0), (7.0/6.0)) * r0 );
    double bb = 1 / ( pow((13.0/7.0), (13.0/6.0)) * r0 );
    
    double Eo = Rep_const / (12*(aa - bb));
    
    double F_max = 50;
    double Mu = 1;
    
    double F_p = 0.15;
    

    
    double phi = 0.0 * M_PI;
    double l_pil = l_min;
    double retraction_period = 10;
    
    normal_distribution<double> T_rev_dist(mean_reversal_period, stdev_reversal_period);
    
    double T_rev_max = mean_reversal_period + 4 * stdev_reversal_period;
    double T_rev_min = max(0.0, mean_reversal_period - 4 * stdev_reversal_period);
    double T_rev;
    
    double P_attach_particle = 0.0;


    double coverage = 0.3;

    int N = round(L*L * coverage)/(l_avg + M_PI * pow(w/2, 2.0));

    double Ls_x = L;
    double Ls_y = L;
    
    double dt_max = 0.1;
    double dt_min = 0.001;
    double dt = dt_max; // this will change to keep the maximum particle movement below dX_max or dTheta_max
    double t = 0.0;
    double tf = 100001;
    cout << tf << endl;
    
    
    double dX_max = 0.1;
    double V_max = 0;
    
    double d_rec_XYTL = 20;
    double d_print_XYTL = 10000;
    double d_check = 1;
    double big_num = (tf + d_rec_XYTL) / dt_min + 1;
    
    double t_last_rec_XYTL = 0;
    double t_last_check = 0;
    double t_last_print_XYTL = 0;
    
    bool check;
    bool rec_XYTL;
    bool print_XYTL = 0;
    int inc = -1;
    
    //save parameter list
    string params = output_directory_lvl_2 + "parameters.txt";
    // writes a text file, saved to input directory, containing the specified string
    ofstream params_txt;
    params_txt.open(params.c_str());
    
    params_txt << "parameter" << "\t" << "value" << "\n"
    << "coverage fraction" << "\t" << coverage << "\n"
    << "number of particles" << "\t" << N << "\n"
    << "mean reversal period" << "\t" << mean_reversal_period << "\n"
    << "standard deviation T_rev" << "\t" << stdev_reversal_period << "\n"
    << "retraction period" << "\t" << retraction_period << "\n"
    << "surface attachment probability" << "\t" << P_attach << "\n"
    << "pilus-particle attachment probability" << "\t" << P_attach_particle << "\n"
    << "rng_seed" << "\t" << seed << "\n"
    << "minimum length" << "\t" << l_min << "\n"
    << "width" << "\t" << w << "\n"
    << "maximum length" << "\t" << l_max << "\n"
    << "side length sorting" << "\t" << s_CL << "\n"
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
    << "list check period" << "\t" << d_check << "\n";
    
    params_txt.close();
    
    
    //initialize output vectors for X, Y, theta
    output_matrix X_Y_T_L(d_print_XYTL, d_rec_XYTL, N, 4);
    
    // initialize bins
    BINS bins(L, L, s_CL);
    
    // initialize N particles
    vector <Particle> particles;
    
    for (int i = 0; i < N; i++)
    {
        double x = ( rng() - 0.5 ) * (Ls_x);
        double y = ( rng() - 0.5 ) * (Ls_y);
        double theta = rng() * 2 * M_PI;
       // double l = l_min + rng() * (l_max - l_min);
        double l = l_avg; // uniform aspect ratio (uniform velocity response)
        
        T_rev = rng() * mean_reversal_period;
        while (T_rev < T_rev_min || T_rev > T_rev_max)
        {T_rev = T_rev_dist(generator);}
        

        Particle p(x, y, l, phi, theta, retraction_period, T_rev, l_pil, F_p, i,
                   N, bins.XI_CL_max, bins.YI_CL_max, L, L);
        
        p.bin(s_CL, L, L);
        
        bins.put_in_bin(p.YI, p.XI, p.ID);
        
        particles.push_back(p);
        
    }
    
    
    
    for (int a = 0; a < mean_reversal_period * 10; a++) //initialize many particles with unsynchronized reversal clocks
    {
        for (int i = 0; i < N; i++)
        {
            
            T_rev = T_rev_dist(generator);
            while (T_rev < T_rev_min || T_rev > T_rev_max)
            {T_rev = T_rev_dist(generator);}

            particles[i].init_reversal_dist(T_rev);
            
        }
    }
    
    
    for (int step = 0; step <= big_num; step++)
    {
        
        check = 0;
        rec_XYTL = 0;
        print_XYTL = 0;
        V_max = 0;
        
        t += dt;
        //cout << step << endl;
        
        if (t > tf){ cout << "tf reached " << t << endl; break; }
        
        if (t - t_last_check > d_check){t_last_check = round(t); check = 1;}
        
        if (t - t_last_rec_XYTL > d_rec_XYTL)
        {
            t_last_rec_XYTL = round(t); rec_XYTL = 1; inc += 1;
            cout << setprecision(9)<< t << "\t" << dt << "\t" << step << "\t" << N << endl;
        }
        
        
        if (t - t_last_print_XYTL > d_print_XYTL)
        {
            t_last_print_XYTL = round(t); print_XYTL = 1;
		cout << "going to print " << t_last_print_XYTL << endl;
        }
        
        
        
        
        for (int i = 0; i < N; i++)
        {
            
            particles[i].move(dt, L, L);
            
            particles[i].reset_calls(N);
            
            particles[i].bin(s_CL, L, L);
            
            bins.put_in_bin(particles[i].YI, particles[i].XI, particles[i].ID);
            
        }
        

        for (int i = 0; i < N; i++)
        {
            // store data for this timestep
            
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
            
            particles[i].Twitch(dt, L, L, s_CL, bins, particles, P_attach, P_attach_particle, T_rev);
            
        }
        
        
        
        bins.clear_bins();
        
        
        for (int i = 0; i < N; i++)
        {
            particles[i].V_tot(Mu);
            V_max = max(V_max, particles[i].V_pole_max);
        }
        
        dt = min(dX_max/V_max, dt_max);
        dt = max(dt, dt_min);
        
        
        if (print_XYTL)
        {
            //write XYTL data to text file, start new matrix for next batch of XYTL data
            string time = ZeroPadNumber((int)round(t), 5);
		cout << time << endl;
            string XYTL_filename = output_directory_lvl_2 + "XYTL_" + "trep_" + time + id_0 + id_1 + ".txt";
            X_Y_T_L.write_to_file(XYTL_filename);

            //initialize output vectors for X, Y, theta
            output_matrix X_Y_T_L(d_print_XYTL, d_rec_XYTL, N, 4);
	    inc = -1;
        }
    }
    

    
    cout << L/2 << endl;
    
    
    
    return 0;
}

