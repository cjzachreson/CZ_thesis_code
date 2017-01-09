//
//  bins_class.h
//  PA_11_01_16
//
//  Created by 11678505 on 11/01/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
// this code was written by Cameron Zachreson (11678505 was my UTS student ID number)
/* This class creates local neighborhoods of particles for efficient location of neighbors 
*/

#ifndef PA_11_01_16_bins_class_h
#define PA_11_01_16_bins_class_h


#endif
// did this change in the original?
class BINS {
    int d1;
    int d2;
    
public:
    vector<vector<vector <int> > > bins;
    int XI_CL_max;
    int YI_CL_max;
    
    BINS(double Lx, double Ly, double s_CL)
    {
        d1 = Ly / s_CL;
        d2 = Lx / s_CL;
        
        bins.resize(d1);
        
        for (int i = 0; i < d1; i++)
        {
            bins[i].resize(d2);
        }
        
        XI_CL_max = Lx / s_CL - 1;
        YI_CL_max = Ly / s_CL - 1;
    }
    
    void put_in_bin(int YI, int XI, int ID)
    {
        bins[YI][XI].push_back(ID);
    }
    
    void clear_bins()
    {
        for (int i = 0; i < d1; i++)
        {
            for (int j = 0; j < d2; j++)
            {
                bins[i][j].resize(0);
            }
        }
    }
    
};
