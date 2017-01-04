//
//  output_matrix_class.h
//  PA_11_01_16
//
//  Created by 11678505 on 11/01/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
//

#ifndef PA_11_01_16_output_matrix_class_h
#define PA_11_01_16_output_matrix_class_h


#endif

class output_matrix{
    int d1;
    int d2;
    int d3;
    
public:
    vector<vector<vector <double> > > XYTL;
    
    //constructor
    output_matrix(double tf, double dt, int N, int num_params)
    {
        d1 = (int)round(tf/dt) + 1;
        d2 = N;
        d3 = num_params;
        this->XYTL.resize(d1);
        
        for (int i = 0; i < d1; i++)
        {
            this->XYTL[i].resize(d2);
            
            for (int j = 0; j < d2; j++){
                
                this->XYTL[i][j].resize(d3);
            }
        }
    }
    
    //updater
    inline void insert_data(int t, int i, int parameter, double val)
    {
        this->XYTL[t][i][parameter] = val;
    }
    
    
    void write_to_file(string filename)
    {
        ofstream XYTL_txt;
        
        XYTL_txt.open(filename.c_str());
        
        cout << "start writing" << endl;
        
        for(int i = 0; i < d2; i++){
            
            for(int t = 0; t < d1; t++)
            {
                XYTL_txt << this->XYTL[t][i][0] << "\t" << this->XYTL[t][i][1]<< "\t" << this->XYTL[t][i][2] << "\t" << this-> XYTL[t][i][3] << "\t";
            }
            
            XYTL_txt << endl;
        }
        XYTL_txt.close();
        cout << "done writing" << endl;
    }
};
