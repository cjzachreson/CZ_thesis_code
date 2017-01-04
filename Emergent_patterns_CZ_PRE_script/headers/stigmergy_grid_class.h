//
//  stigmergy_gird_class.h
//  PA_11_01_16
//
//  Created by 11678505 on 11/01/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
//

#ifndef PA_11_01_16_stigmergy_gird_class_h
#define PA_11_01_16_stigmergy_gird_class_h


#endif

using namespace std;

class Stigmergy_Grid{
    int d1;
    int d2;
    double P_min;
    double P_max;
    vector <int> linear_indices; // will be filled with the indices for all particles
    vector <int> seg_ind_x; // indices for a rod segment on the origin
    vector <int> seg_ind_y;
    double diff_x_s;
    double diff_y_s;
    
    double diff_x_l;
    double diff_y_l;
    
    double gamma_s;
    double gamma_l;
    
    double F_seg_x;
    double F_seg_y;
    
    
    
    void wrap_pixel(int & row, int & column, int xi_max, int yi_max) const
    {
        if (row > yi_max) {row -= yi_max + 1;}
        if (row < 0) {row += yi_max + 1;}
        if (column > xi_max) {column -= xi_max + 1;}
        if (column < 0) {column += xi_max + 1;}
    }
    
public:
    vector<vector <double> > Stig_Grid_P_now;
    vector<vector <double> > Stig_Grid_P_next;
    
    vector<vector <double> > Stig_Grid_S_now;
    vector<vector <double> > Stig_Grid_S_next;
    
    vector<vector <double> > Stig_Grid_L_now;
    vector<vector <double> > Stig_Grid_L_next;
    
    int XI_SG_max;
    int YI_SG_max;
    
    double F_env_x;
    double F_env_y;
    double torque_env;
    
    
    
    Stigmergy_Grid(double Lx, double Ly, double s_SG, double w, double p_min, double p_max, double gam_s, double gam_l)
    {
        d1 = round(Ly / s_SG);
        d2 = round(Lx / s_SG);
        P_min = p_min;
        P_max = p_max;
        gamma_s = gam_s;
        gamma_l = gam_l;
        
        this->Stig_Grid_P_now.resize(d1);
        this->Stig_Grid_S_now.resize(d1);
        this->Stig_Grid_L_now.resize(d1);
        
        for (int i = 0; i < d1; i++)
        {
            this->Stig_Grid_P_now[i].resize(d2);
            this->Stig_Grid_S_now[i].resize(d2);
            this->Stig_Grid_L_now[i].resize(d2);
            
        }
        
        this->Stig_Grid_P_next = this->Stig_Grid_P_now;
        this->Stig_Grid_S_next = this->Stig_Grid_S_now;
        this->Stig_Grid_L_next = this->Stig_Grid_L_now;
        
        XI_SG_max = round(Lx / s_SG) - 1;
        YI_SG_max = round(Ly / s_SG) - 1;
        
        int rad = (int)floor(((w*1.5)/2.0) / s_SG);
        vector <int> seg_diam;
        vector <int> Y_seg_top_edge;
        vector <int> Y_seg_bottom_edge;
        
        for (int i = -rad; i<(rad+1); i++)
        {
            seg_diam.push_back(i);
            
            Y_seg_top_edge.push_back(floor(sqrt(rad*rad - seg_diam.back()*seg_diam.back())));
            
            Y_seg_bottom_edge.push_back(- Y_seg_top_edge.back());
            
            for (int j = Y_seg_bottom_edge.back(); j <= Y_seg_top_edge.back(); j++ )
            {
                seg_ind_y.push_back(j);
                seg_ind_x.push_back(i);
            }
        }
    }
    
    void update_stig_grid()
    {
        this->Stig_Grid_P_now = this->Stig_Grid_P_next;
        this->Stig_Grid_S_now = this->Stig_Grid_S_next;
        this->Stig_Grid_L_now = this->Stig_Grid_L_next;
        
    }
    
    void clear_for_next_particle()
    {
        this->torque_env = 0;
        this->F_env_x = 0;
        this->F_env_y = 0;
        
    }
    
    void index_particle_to_grid(int x_i, int y_i, double theta_i, double l_i, double w, double s_SG, double Lx, double Ly)
    {
        // x and y coords for each segment as if the rod was on the origin with at angle of 0 degrees
        double l_seg = floor(l_i/(w * 0.5));
        
        double centroid_X;
        double centroid_Y;
        
        for (int i = 0; i<=(l_seg); i++)
        {
            centroid_X = -l_i/2 + l_i/l_seg * i;
            centroid_Y = 0;
            
            double r_seg = centroid_X;
            
            rotate_coords_2D(centroid_X, centroid_Y, theta_i);
            
            centroid_X = (int)round(centroid_X/s_SG);
            centroid_Y = (int)-round(centroid_Y/s_SG);
	    
            diff_x_s = 0;
            diff_y_s = 0;
            
            diff_x_l = 0;
            diff_y_l = 0;
            
            F_seg_x = 0;
            F_seg_y = 0;
            
            
            for (int j = 0; j<(int)seg_ind_y.size(); j++)
            {
                int row = seg_ind_y[j] + centroid_Y + y_i;
                int column = seg_ind_x[j] + centroid_X + x_i;
		
               
		
                wrap_pixel(row, column, XI_SG_max, YI_SG_max);
                
                int row_p1 = row + 1;
                int row_m1 = row - 1;
                
                int column_p1 = column + 1;
                int column_m1 = column - 1;
                
                wrap_pixel(row_p1, column_p1, XI_SG_max, YI_SG_max);
                wrap_pixel(row_m1, column_m1, XI_SG_max, YI_SG_max);
                
                diff_x_s += (Stig_Grid_S_now[row][column_m1] - Stig_Grid_S_now[row][column_p1]) / 2;
                
                diff_y_s += (Stig_Grid_S_now[row_p1][column] - Stig_Grid_S_now[row_m1][column]) / 2;
                
                diff_x_l += (Stig_Grid_L_now[row][column_m1] - Stig_Grid_L_now[row][column_p1]) / 2;
                
                diff_y_l += (Stig_Grid_L_now[row_p1][column] - Stig_Grid_L_now[row_m1][column]) / 2;
                
                linear_indices.push_back(x_y_to_lind(d1, column, row));
            }
            
            F_seg_x = -gamma_s * diff_x_s - gamma_l * diff_x_l;
            F_seg_y = -gamma_s * diff_y_s - gamma_l * diff_y_l;
            
            double F_seg_mag = sqrt(F_seg_x * F_seg_x + F_seg_y  * F_seg_y);
            
            
            if ((double)round(F_seg_mag*10000)/10000 != 0)
            {
                double U_F_seg_x = F_seg_x / F_seg_mag;
                double U_F_seg_y = F_seg_y / F_seg_mag;
                
                double U_perp_seg = cos(theta_i) * U_F_seg_y - sin(theta_i) * U_F_seg_x;
                
                double torque_seg = r_seg * F_seg_mag * U_perp_seg;
                
                this->torque_env += torque_seg;
                
                this -> F_env_x += F_seg_x;
                this -> F_env_y += F_seg_y;
                
            }
        }
    }
    
    double find_P_attch(int row, int column, double s_SG) const
    {
        wrap_pixel(row, column, XI_SG_max, YI_SG_max);
        double p_attch = (P_max - P_min) * (this->Stig_Grid_P_now[row][column] / (s_SG * s_SG)) + P_min;
        return p_attch;
    }
    
    void decay(double beta_p, double beta_s, double beta_l, double dt)
    {
        for (int i = 0; i < d1; i++)
        {
            for (int j = 0; j < d2; j++)
            { this->Stig_Grid_P_next[i][j] -= this->Stig_Grid_P_next[i][j] * beta_p * dt;
                this->Stig_Grid_S_next[i][j] -= this->Stig_Grid_S_next[i][j] * beta_s * dt;
                this->Stig_Grid_L_next[i][j] -= this->Stig_Grid_L_next[i][j] * beta_l * dt;
              
            }
        }
    }
    
    void add_counts(double k_p, double k_s, double k_l, double dt, double c_max)
    {
        int d = (int) linear_indices.size();
        
        sort(linear_indices.begin(), linear_indices.end());
        
        linear_indices.push_back(-1);
        
        for (int i = 0; i < d; i++)
        {   if (linear_indices[i] != linear_indices[i + 1])
        {
            int column = x_from_lind(d1, linear_indices[i]);
            int row = y_from_lind_and_x_ind(d1, column, linear_indices[i]);
            
            this->Stig_Grid_P_next[row][column] +=
            k_p * ((c_max - this->Stig_Grid_P_now[row][column]) / c_max) * dt;
            
            this->Stig_Grid_S_next[row][column] +=
            k_s * ((c_max - this->Stig_Grid_S_now[row][column]) / c_max) * dt;
            
            this->Stig_Grid_L_next[row][column] +=
            k_l * ((c_max - this->Stig_Grid_L_now[row][column]) / c_max) * dt;
        }
        }
        
        linear_indices.resize(0);
    }
};
