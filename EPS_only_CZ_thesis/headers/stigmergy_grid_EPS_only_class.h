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
    
    
    int XI_SG_max;
    int YI_SG_max;
    
    
    
    
    Stigmergy_Grid(double Lx, double Ly, double s_SG, double w, double p_min, double p_max)
    {
        d1 = round(Ly / s_SG);
        d2 = round(Lx / s_SG);
        P_min = p_min;
        P_max = p_max;
        
        this->Stig_Grid_P_now.resize(d1);

        for (int i = 0; i < d1; i++)
        {
            this->Stig_Grid_P_now[i].resize(d2);
            
        }
        
        this->Stig_Grid_P_next = this->Stig_Grid_P_now;
        
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

            for (int j = 0; j<(int)seg_ind_y.size(); j++)
            {
                int row = seg_ind_y[j] + centroid_Y + y_i;
                int column = seg_ind_x[j] + centroid_X + x_i;

                wrap_pixel(row, column, XI_SG_max, YI_SG_max);
                
                linear_indices.push_back(x_y_to_lind(d1, column, row));
            }
            
        }
    }
    
    double find_P_attch(int row, int column, double s_SG) const
    {
        wrap_pixel(row, column, XI_SG_max, YI_SG_max);
        double p_attch = (P_max - P_min) * (this->Stig_Grid_P_now[row][column] / (s_SG * s_SG)) + P_min;
        return p_attch;
    }
    
    void decay(double beta_p, double dt)
    {
        for (int i = 0; i < d1; i++)
        {
            for (int j = 0; j < d2; j++)
            { this->Stig_Grid_P_next[i][j] -= this->Stig_Grid_P_next[i][j] * beta_p * dt;

              
            }
        }
    }
    
    void add_counts(double k_p, double dt, double c_max)
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
            
        }
        }
        
        linear_indices.resize(0);
    }
};
