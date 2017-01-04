//
//  image_matrix_class.h
//  PA_11_01_16
//
//  Created by 11678505 on 11/01/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
//

#ifndef PA_11_01_16_image_matrix_class_h
#define PA_11_01_16_image_matrix_class_h


#endif

// only works if pngwriter library is installed in this directory, otherwise use text image
//#include "/usr/local/include/pngwriter.h"

class image_matrix{
    
public:
    
    
    //constructor
    
    image_matrix(string filename, vector<vector<double> > & grid, int XI_max, int YI_max)
    {
        ofstream image_file;
        image_file.open(filename.c_str());

        for (int i = 0; i < YI_max; i++)
        {
            for (int j = 0; j < XI_max; j++)
            {
                image_file << grid[i][j] << "\t";
                
            }
            
            image_file << endl;
        }
        
    }
    
    
    /*   image_matrix(string filename, vector<vector<double> > & grid, int XI_max, int YI_max)
     {
     pngwriter png(YI_max, XI_max, 1, filename.c_str());
     
     for (int i = 0; i < YI_max; i++)
     {
     for (int j = 0; j < XI_max; j++)
     {
     png.plot(j, YI_max - i, 1 - grid[i][j], 1 - grid[i][j], 1 - grid[i][j]);
     }
     }
     png.close();
     }
     */
    
    
    
    
};
