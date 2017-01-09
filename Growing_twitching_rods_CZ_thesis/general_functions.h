//
//  general_functions.h
//  PA_11_01_16
//
//  Created by 11678505 on 11/01/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
// this code was written by Cameron Zachreson (11678505 was my UTS student ID number)
/* This class contains general use functions
*/

#ifndef PA_11_01_16_general_functions_h
#define PA_11_01_16_general_functions_h


#endif

#include <random>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

inline double rng()
{
    /*
     // Christian Wolff: Adopted from the old implementation of rand(3) in glibc.
     // By no means acceptable for cryptography, but good enough for our
     // application and fast.
     this->rng_register_ = ((this->rng_register_ * 1103515245) + 12345) & 0x7fffffff;
     return ((double) this->rng_register_) / ((double) 0x80000000);
     */
    return ((double) rand() / (double) RAND_MAX);
}

int index_X_to_grid(double side_length, double Lx, double X_coord)

{
    int X_index = floor(X_coord / side_length) + (Lx / side_length / 2);
    return X_index;
}


int index_Y_to_grid(double side_length, double Ly, double Y_coord)
{
    int Y_index = abs(ceil(Y_coord / side_length) - (Ly / side_length / 2));
    return Y_index;
}

int x_y_to_lind(int d1, int column, int row)
{
    int lind = (column + 1) * d1 - (d1 - row);
    return lind;
}

int x_from_lind(int d1, int lind)
{
    int x_ind = (int)floor(lind/d1);
    return x_ind;
}

int y_from_lind_and_x_ind(int d1, int x_ind, int lind)
{
    int y_ind = lind - x_ind * d1;
    return y_ind;
}

int sign(double x)
{
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}


double dist_2D_PB(double X1, double X2, double Y1, double Y2, double Lx, double Ly)
{
    double d1 = min(abs(X1-X2), Lx - abs(X1 - X2));
    double d2 = min(abs(Y1-Y2), Ly - abs(Y1 - Y2));
    double dist = sqrt((d1)*(d1) + (d2)*(d2));
    return dist;
}

double dist_2D(double dX, double dY)
{
    double dist = sqrt((dX)*(dX) + (dY)*(dY));
    return dist;
}


double wrap_point(double p, double L)
{
    if (abs(p) > L/2)
    {double p1 = (L - abs(p)) * -sign(p);
        return p1;}
    else return p;
}



double wrap_vec(double p1, double p2, double L)
{
    double dist = p2 - p1;
    double dist_P = min(abs(dist), L - abs(dist)) * sign(dist);
    if (abs(dist_P) < abs(dist)) { dist_P = -dist_P; }
    return dist_P;
}


void rotate_coords_2D(double & coord_x, double & coord_y, const double theta)
{
    double x_rot = coord_x * cos(theta) + coord_y * -sin(theta);;
    coord_y = coord_x * sin(theta) + coord_y * cos(theta);
    coord_x = x_rot;
}


string ZeroPadNumber(int num, int digits)
{
    std::ostringstream ss;
    ss << std::setw( digits ) << std::setfill( '0' ) << num;
    return ss.str();
}

string num2str(double num)
{
    ostringstream ss;
    ss << num;
    return ss.str();
}

