//
//  main.cpp
//  growing twitching rods
//
//  Created by 11678505 on 21/04/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
//
// this code was written by Cameron Zachreson (11678505 was my UTS student ID number)
/* This script simulates surface-motile bacteria that grow (at a constant rate) and divide upon reaching a critical length
the code outputs a matrix of position, orientation, and length data for N bacteria in the following format:
X1(t = 1), Y1(t = 1); theta1(t = 1); length1(t = 1); X1(t = 2); Y1(t = 2); theta1(t = 2); length1(t = 2);... length1(tf)
X2(t = 1), Y2(t = 1); theta2(t = 1); length2(t = 1); X2(t = 2); Y2(t = 2); theta2(t = 2); length2(t = 2);... length2(tf)
.
.
.
XN(t = 1), YN(t = 1); thetaN(t = 1); lengthN(t = 1); XN(t = 2); YN(t = 2); thetaN(t = 2); lengthN(t = 2);... lengthN(tf)
this matrix can be read by the XYTL class for analysis in the .m files provided.
the script also creates text files containing the values of each pixel in the stigmergy grid 
the frequency of data recording can be set by changing d_rec
*/
#include <cstdlib>
#include <algorithm>
#include <random>

#include <iostream>
#include <sys/stat.h>
#include <sstream>

#include <array>



#include "general_functions.h"
#include "bins_class.h"
#include "particle_class.h"
#include "output_matrix_class.h"


int main(int argc, const char * argv[]) {
    
    double seed = 1;
    
    string state = "growing_rods_tst_2";
    
    double growth_rate = 0.001;
    
    double mean_reversal_period = 1000;//need to tune this for WT behavior
    
    double P = 0.1; // attachment probability, set P = 0 for nonmotile case

    
    
    string output_directory = "/Users/11678505/Desktop/growing_twitching_rods_output/" + state;
    mkdir(output_directory.c_str(), 0777);
    output_directory = output_directory + "/output/";
    mkdir(output_directory.c_str(), 0777);

    default_random_engine generator;
    generator.seed(seed); // for gaussian distribution (reversal period)
    srand(seed);//for uniform distributions
    
    
    ostringstream label_lvl_1;
    label_lvl_1 << "_tr_" << mean_reversal_period;
    
    ostringstream label_lvl_2;
    label_lvl_2 << "_GR_" << growth_rate;
    
    ostringstream label_lvl_3;
    label_lvl_3 << "_Pa_" << P;
    
    string id_1 = label_lvl_1.str();
    string id_2 = label_lvl_2.str();
    string id_3 = label_lvl_3.str();
    
    string id = id_1 + id_2 + id_3;
    
    string output_directory_lvl_2 = output_directory + state + id + "/";
    mkdir(output_directory_lvl_2.c_str(), 0777);
    

    double r0 = 1;
    
    double Rep_const = 1;
    
    double aa = 1 / ( pow((13.0/7.0), (7.0/6.0)) * r0 );
    double bb = 1 / ( pow((13.0/7.0), (13.0/6.0)) * r0 );
    
    double Eo = Rep_const / (12*(aa - bb));
    
    double F_max = 20;
    double Mu = 1;
    
    double F_p = 1.5;
    
    double dl = 0.0;
    
    double l_min = 3.0;
    double w = 1;
    double l_max = 2 * l_min + w;
    double s_CL = l_max + w;
    
    double L = s_CL * 200;
    double Ls_x = 0;
    double Ls_y = 0;
    
    double X_max = 0;
    double Y_max = 0;
    
    double phi = 0.5 * M_PI;
    double l_pil = l_min;
    double retraction_period = 10;
    
    //normal_distribution<double> T_rev_dist(mean_reversal_period, stdev_reversal_period);
    // this version of the code used a uniform distribution of reversal frequencies
    double T_rev_max = mean_reversal_period * 2; //+ 4 * stdev_reversal_period; 
    double T_rev_min = 0.0;
    double T_rev;
    
    double P_attach_particle = 0; // pilus cross-linking is disabled
    int N = 1;
    int N_new = 0;
    
    double pop_rn;
    int sub_pop;
    
    double dt_max = 0.1;
    double dt_min = 0.005;
    double dt = dt_max; // this will change to keep the maximum particle movement below dX_max or dTheta_max
    double t = 0.0;
    double replication_cycles = 14;
    double tf = round((l_max - l_min) / growth_rate) * replication_cycles;
    cout << tf << endl;
    double big_num = tf / dt_min + 1;
    
    double dX_max = 0.1;
    double V_max = 0;
    
    double d_rec_XYTL = 100;
    double d_check = 1;
    
    double t_last_rec_XYTL = 0;
    double t_last_check = 0;
    
    bool check;
    bool rec_XYTL;
    int inc = -1;
    
    //save parameter list
    string params = output_directory_lvl_2 + "parameters.txt";
    // writes a text file, saved to input directory, containing the specified string
    ofstream params_txt;
    params_txt.open(params.c_str());
    
    params_txt << "parameter" << "\t" << "value" << "\n"
    << "number of particles" << "\t" << N << "\n"
    << "mean reversal period" << "\t" << mean_reversal_period << "\n"
    << "standard deviation T_rev" << "\t" << "uniform" << "\n"
    << "retraction period" << "\t" << retraction_period << "\n"
    << "attachment probability" << "\t" << P << "\n"
    << "pilus-particle attachment probability" << "\t" << P_attach_particle << "\n"
    << "rng_seed" << "\t" << seed << "\n"
    << "growth_rate" << "\t" << growth_rate << "\n"
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
    << "max possible time step" << "\t" << dt_max << "\n"
    << "min possible time step" << "\t" << dt_min << "\n"
    << "t_final" << "\t" << tf << "\n"
    << "dX max" << "\t" << dX_max << "\n"
    << "config recording period" << "\t" << d_rec_XYTL << "\n"
    << "list check period" << "\t" << d_check << "\n";
    
    params_txt.close();
    
    string XYTL_filename = output_directory_lvl_2 + "XYTL" + id + ".txt";
    
    
    //initialize output vectors for X, Y, theta
    
    output_matrix X_Y_T_L(tf, d_rec_XYTL, N, 4);
    
    // initialize bins
    BINS bins(L, L, s_CL);
    
    // initialize N particles
    vector <Particle> particles;
    
    for (int i = 0; i < N; i++)
    {
        double x = ( rng() - 0.5 ) * (Ls_x);
        double y = ( rng() - 0.5 ) * (Ls_y);
        double theta = rng() * 2 * M_PI;
        double l = l_min + rng() * (l_max - l_min);
        
        
        T_rev = rng() * T_rev_max; // uniform distribution
        //while (T_rev < T_rev_min || T_rev > T_rev_max) // gaussian distribution
        //{T_rev = T_rev_dist(generator);}
        
        Particle p(x, y, l, phi, theta, retraction_period, T_rev, l_pil, F_p, i,
                   N, bins.XI_CL_max, bins.YI_CL_max, L, L, growth_rate);
        
        p.bin(s_CL, L, L);
        
        bins.put_in_bin(p.YI, p.XI, p.ID);
        
        particles.push_back(p);
        
    }
    
    
    
    /* for (int a = 0; a < mean_reversal_period * 10; a++) //initialize many particles with unsynchronized reversal clocks
     {
     for (int i = 0; i < N; i++)
     {
     
     T_rev = T_rev_dist(generator);
     while (T_rev < T_rev_min || T_rev > T_rev_max)
     {T_rev = T_rev_dist(generator);}
     
     particles[i].init_reversal_dist(T_rev);
     
     }
     }*/
    
    
    for (int step = 0; step < big_num; step++)
    {
        
        check = 0;
        rec_XYTL = 0;
        V_max = 0;
        N_new = 0;
        
        t += dt;
        //cout << step << endl;
        
        if (t >= tf){break;}
        if (abs(Y_max) > (L/2 - s_CL)) {break;}
        if (abs(X_max) > (L/2 - s_CL)) {break;}
        
        if (t - t_last_check > d_check){t_last_check = round(t); check = 1;}
        
        if (t - t_last_rec_XYTL > d_rec_XYTL)
        {
            t_last_rec_XYTL = round(t); rec_XYTL = 1; inc += 1;
            cout << setprecision(9)<< t << "\t" << dt << "\t" << step << "\t" << N << endl;
        }
        
        
        
        for (int i = 0; i < N; i++)
        {
            
            particles[i].move(dt, L, L);
            
            if (abs(particles[i].X) > X_max){X_max = abs(particles[i].X);}
            if (abs(particles[i].Y) > Y_max){Y_max = abs(particles[i].Y);}
            
            particles[i].grow(dt);
            
            if (particles[i].l > l_max)
            {
                
                N_new += 1;
                
                int id = N + N_new - 1;
                
                dl = rng();
                
                double x = particles[i].X - (particles[i].l/4 + w/2) * cos(particles[i].theta);
                double y = particles[i].Y - (particles[i].l/4 + w/2) * sin(particles[i].theta);
                double theta = particles[i].theta;
                double l = (particles[i].l - w)/2 - dl;
                
                
                T_rev = rng() * T_rev_max;
                //while (T_rev < T_rev_min || T_rev > T_rev_max)
                //{T_rev = T_rev_dist(generator);}
                
                
                Particle p(x, y, l, phi, theta, retraction_period, T_rev, l_pil, F_p, id,
                           N, bins.XI_CL_max, bins.YI_CL_max, L, L, growth_rate);
                
                particles.push_back(p);
                
                particles[i].divide(L, L, dl); //alters position and length of mother particle
                
            }
            
        }
        
        
        N += N_new;
        if (N_new > 0)
        {
            X_Y_T_L.more_particles(N);
            check = 1;
        }
        
        for (int i = 0; i < N; i++)
        {
            particles[i].reset_calls(N);
        }
        

        for (int i = 0; i < N; i++)
        {
            
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
            
            T_rev = rng() * T_rev_max;//T_rev_dist(generator);
            /*while (T_rev < T_rev_min || T_rev > T_rev_max)
             {
             T_rev = T_rev_dist(generator);
             }*/
            
            particles[i].Twitch(dt, L, L, s_CL, bins, particles, P, P_attach_particle, T_rev);
            
        }
        
        
        
        bins.clear_bins();
        
        
        for (int i = 0; i < N; i++)
        {
            particles[i].V_tot(Mu);
            V_max = max(V_max, particles[i].V_pole_max);
        }
        
        dt = min(dX_max/V_max, dt_max);
        dt = max(dt, dt_min);
        
        
    }
    
    //write data to text file
    
    X_Y_T_L.write_to_file(XYTL_filename);
    
    cout << L/2 << endl;
    
    
    
    return 0;
}

