//
//  particle_class.h
//  PA_11_01_16
//
//  Created by 11678505 on 11/01/2016.
//  Copyright (c) 2016 11678505. All rights reserved.
//

#ifndef PA_11_01_16_particle_class_h
#define PA_11_01_16_particle_class_h


#endif

class Particle {
    double w;
    double vx;
    double vy; // transltional velocity of this particle
    double v_Ang; // angular velocity of this rod shaped particle
    double phi; // angular range of pilus
    double pilus_range;//the radius of attachment for a particular attachment
    double pilus_length; // the maximum possible attachment radius
    double attachment_point_x;
    double attachment_point_y;// x and y coordinates for the attachment point of this one's pilus
    double retraction_timer;
    double retraction_period;
    double reversal_timer;
    double F_p; // retraction force
    bool attachment_surface; // pilus attached to the substratum?
    bool attachment_particle; // pilus attached to another particle?
    vector <int> hood; // particles near this one
    vector <int> hood_ap; //particles near this one's attachment point
    int XI_max; // maximum possible XY index values in sorting lattice
    int YI_max;
    int j_ap; // the id of the particle i attaches to with its pilus
    double ap_o_par_j; // the lever arm distance on particle j attached via particle i's pilus
    double ap_o_perp_j; // the perpendicular component of the attachment point relative to particle j's midline
    double GR;
    
    void find_hood_generic(const BINS &, vector <int> &, int, int);
    void clear_hood(){hood.resize(0);}
    void clear_hood_ap(){hood_ap.resize(0);}
    
    
public:
    
    double X;
    double Y;
    double l;
    int XI;
    int YI;
    int XI_ap;
    int YI_ap;
    int xi_sg;
    int yi_sg;
    int ID;
    double theta;
    double F_x;
    double F_y;
    double torque;
    double V_pole_max;
    double P_attch;
    vector <bool> calls;
    
    // put reversal distribution inside the particle class as a characteristic, then initialization of a particle can just pick between the different distributions for establishing T_rev
    
    Particle(double, double, double, double, double, double, double, double, double, int, int, int, int, double, double, double);
    
    void init_reversal_dist(double);
    
    void find_hood(int, int, const BINS &);
    
    void Repulsion(vector<Particle> &,  const BINS &, int, double, double, double, double, double, bool);
    
    void Twitch(const double, const double, const double, const double, const BINS &, vector<Particle> &, const double, const double, double);
    
    void move(double dt, double Lx, double Ly)
    {
        this->X += vx * dt;
        this->Y += vy * dt;
        this->theta += v_Ang * dt;
        
        this->X = wrap_point(this->X, Lx);
        this->Y = wrap_point(this->Y, Ly);
    }
    
    void reset_calls(int N)
    {
        this -> calls.assign(N, true);
    }

    
    void grow(double dt)
    {
        this -> l += GR * dt;
    }
    
    void divide(double Lx, double Ly, double dl)
    {
        
        this->X = this->X + (this->l/4 + w/2) * cos(this->theta);
        this->Y = this->Y + (this->l/4 + w/2) * sin(this->theta);
        this->l = (this->l - w)/2 + dl;
        
        this->X = wrap_point(this->X, Lx);
        this->Y = wrap_point(this->Y, Ly);
        
    }
    
    
    void V_tot(double Mu)
    {
        
        vx = this->F_x / (Mu * this->l);
        vy = this->F_y / (Mu * this->l);
        v_Ang =(12 * this->torque)/(Mu * pow(this->l, 3.0));

        this->V_pole_max = dist_2D(vx, vy) + abs((this->l/2) * v_Ang);
        
        this->F_x = 0.0;
        this->F_y = 0.0;
        this->torque = 0.0;
 
    }
    
    
    void bin(double s_CL, double Lx, double Ly)
    {
        this->XI = index_X_to_grid(s_CL, Lx, this->X);
        this->YI = index_Y_to_grid(s_CL, Ly, this->Y);
    }
    
    void bin_ap(double s_CL, double Lx, double Ly, double ap_x, double ap_y)
    {
        this->XI_ap = index_X_to_grid(s_CL, Lx, ap_x);
        this->YI_ap = index_Y_to_grid(s_CL, Ly, ap_y);
    }
    
    void anneal(double Lx, double Ly, double N, double dt)
    {
        double theta_i = this -> theta;
        double X_i = this->X;
        double Y_i = this->Y;
        
        this->X = X_i + (rng() - 0.5)*2 * dt;
        this->Y = Y_i + (rng() - 0.5)*2 * dt;
        this->theta = theta_i + (rng() - 0.5) * dt;
        
        this->X = wrap_point(this->X, Lx);
        this->Y = wrap_point(this->Y, Ly);
        
        this->F_x = 0.0;
        this->F_y = 0.0;
        this->torque = 0.0;
        this->calls.assign(N, true);
    }
};

Particle::Particle(double x_i, double y_i, double l_i, double phi_i, double theta_i,
                   double T_ret, double T_rev, double l_pil, double pilus_force, int I,
                   int N, int xi_max, int yi_max, double Lx, double Ly, double growth_rate)
{
    

    
   
    
    v_Ang = 0.0;
    XI_max = xi_max;
    YI_max = yi_max;
    this->X = x_i;
    this->Y = y_i;
    this->X = wrap_point(this->X, Lx);
    this->Y = wrap_point(this->Y, Ly);
    
    this->l = l_i;
    w = 1;
    this->vx = 0.0;
    this->vy = 0.0;
    this->ID = I;
    this->theta = theta_i;
    this->calls.assign(N, true);
    this->F_x = 0.0;
    this->F_y = 0.0;
    this->torque = 0.0;
    phi = phi_i;
    pilus_range = l_pil;
    retraction_period = T_ret;
    retraction_timer = rng() * T_ret;
    reversal_timer = T_rev;
    F_p = pilus_force;
    attachment_surface = 0;
    attachment_particle = 0;
    GR = growth_rate;
    
    //cout << this-> ID << "\t"<< "created, torque =" << "\t" << this-> torque << endl;
}

void Particle::init_reversal_dist(double t_rev)
{
    this->reversal_timer -= 1;
    
    if (this->reversal_timer < 0)
    {this->reversal_timer = t_rev;}
}


void Particle::Twitch(const double dt, const double Lx, const double Ly, const double s_CL, const BINS & bins, vector<Particle> & p, const double P_attach_surface, const double P_attch_particle, double T_rev)
{
    retraction_timer -= dt;
    reversal_timer -= dt;
    
    if (retraction_timer <= 0)
    {
        retraction_timer = rng() * retraction_period;
        
        if (reversal_timer <= 0 )
        {this->theta += M_PI; reversal_timer = T_rev;}
        
        double pilus_angle = (rng() - 0.5) * phi;
        pilus_length = rng() * pilus_range;
        
        attachment_point_x = (this->l + w) / 2 + pilus_length * cos(pilus_angle);
        attachment_point_y = pilus_length * sin(pilus_angle);
        
        rotate_coords_2D(attachment_point_x, attachment_point_y, this->theta);
        
        attachment_point_x += this->X;
        attachment_point_y += this->Y;
        
        attachment_point_x = wrap_point(attachment_point_x, Lx);
        attachment_point_y = wrap_point(attachment_point_y, Ly);
        
        
        // surface attachment or particle-particle attachment?
        // first find the particles near the attachment point - ie the neighborhood of the attachment point
        bin_ap(s_CL, Lx, Ly, attachment_point_x, attachment_point_y);
        
        clear_hood_ap();
        //find_hood_attachment(bins, hood_ap);
        find_hood_generic(bins, hood_ap, XI_ap, YI_ap);
        
        // iterate through to see if attachment point is on the body of another particle
        
        int N_hood_ap = (int)hood_ap.size();
        
        for (int i = 0; i < N_hood_ap; i++)
        {
            attachment_particle = 0;
            
            int j = hood_ap[i];
            
            double dX_api_cj = -wrap_vec(attachment_point_x, p[j].X, Lx);
            double dY_api_cj = -wrap_vec(attachment_point_y, p[j].Y, Ly);
            
            double d_api_cj = dist_2D(dX_api_cj, dY_api_cj);
            
            double Ux_api_cj = dX_api_cj / d_api_cj;
            double Uy_api_cj = dY_api_cj / d_api_cj;
            
            if (d_api_cj <= (p[j].l + w)/2)
            {
                double d_api_backbone_j =
                d_api_cj * (Uy_api_cj * cos(p[j].theta) - Ux_api_cj * sin(p[j].theta));
                
                if (abs(d_api_backbone_j) <= (w/2))
                {
                    double rlev_api_j =
                    d_api_cj * (Ux_api_cj * cos(p[j].theta) + Uy_api_cj * sin(p[j].theta));
                    
                    if (abs(rlev_api_j) > (p[j].l/2) && dist_2D(d_api_backbone_j, abs(rlev_api_j)-p[j].l/2) > (w/2))
                    {
                        attachment_particle = 0;
                    }
                    
                    else
                    {
                        attachment_particle = rng() > (1 - P_attch_particle);
                        
                        if (attachment_particle)
                        {
                            j_ap = p[j].ID; // the id of the particle that i attaches to with its pilus
                            ap_o_par_j = rlev_api_j;
                            ap_o_perp_j = d_api_backbone_j;
                            
                            break;
                        }
                    }
                    
                }
            }
            
            
            
        }
        
        
        if (attachment_particle == 0)
        {
            this->P_attch = P_attach_surface;
            
            attachment_surface = rng() > (1 - P_attch);
        }
        
    }
    
    if (attachment_particle)
    {
        double pole_x = ((this->l + w) / 2) * cos(this->theta) + this->X;
        double pole_y = ((this->l + w) / 2) * sin(this->theta) + this->Y;
        
        pole_x = wrap_point(pole_x, Lx);
        pole_y = wrap_point(pole_y, Ly);
        
        int j = j_ap;
        attachment_point_x = ap_o_par_j;
        attachment_point_y = ap_o_perp_j;
        
        rotate_coords_2D(attachment_point_x, attachment_point_y, p[j].theta);
        
        attachment_point_x = wrap_point(attachment_point_x + p[j].X, Lx);
        attachment_point_y = wrap_point(attachment_point_y + p[j].Y, Ly);
        
        double d_pil = dist_2D_PB(pole_x, attachment_point_x, pole_y, attachment_point_y, Lx, Ly);
        
        double twitch_vec_x = -wrap_vec(attachment_point_x, pole_x, Lx) / d_pil;
        double twitch_vec_y = -wrap_vec(attachment_point_y, pole_y, Ly) / d_pil;
        
        double twitch_vec_perp_i = (cos(this->theta) * twitch_vec_y) - (sin(this->theta) * twitch_vec_x);
        
        this -> F_x += F_p/2 * twitch_vec_x;
        this -> F_y += F_p/2 * twitch_vec_y;
        this -> torque += F_p/2 * twitch_vec_perp_i * (this->l + w)/2;
        
        double twitch_vec_perp_j =   (cos(p[j].theta) * -twitch_vec_y) - (sin(p[j].theta) * -twitch_vec_x);
        
        p[j].F_x -= F_p/2 * twitch_vec_x;
        p[j].F_y -= F_p/2 * twitch_vec_y;
        p[j].torque += F_p/2 * twitch_vec_perp_j * ap_o_par_j;
        
    }
    
    
    if (attachment_surface)
    {
        double pole_x = ((this->l + w) / 2) * cos(this->theta) + this->X;
        double pole_y = ((this->l + w) / 2) * sin(this->theta) + this->Y;
        
        pole_x = wrap_point(pole_x, Lx);
        pole_y = wrap_point(pole_y, Ly);
        
        double d_pil = dist_2D_PB(pole_x, attachment_point_x, pole_y, attachment_point_y, Lx, Ly);
        
        double twitch_vec_x = -wrap_vec(attachment_point_x, pole_x, Lx) / d_pil;
        double twitch_vec_y = -wrap_vec(attachment_point_y, pole_y, Ly) / d_pil;
        
        double twitch_vec_perp = cos(this->theta) * twitch_vec_y - sin(this->theta) * twitch_vec_x;
        
        this -> F_x += F_p * twitch_vec_x;
        this -> F_y += F_p * twitch_vec_y;
        this -> torque += F_p * twitch_vec_perp * (this->l + w)/2;
    }
    
   // cout << this-> ID << "\t"<< "twitched, torque =" << "\t" << this-> torque << endl;
}


void Particle::Repulsion(vector<Particle> & particles, const BINS & b, int i, double r0, double Eo, double F_max, double Lx, double Ly, bool hood_check)
{
    if (hood_check)
    {
        clear_hood();
        
        find_hood_generic(b, hood, XI, YI);
    }
    
    int N_hood = (int)hood.size();
    vector <int> new_hood;
    
    
    for (int jI = 0; jI < N_hood; jI++)
    {
        int j = hood[jI];
        bool call_ij = this->calls[j] & particles[j].calls[i];
        
        if (j != i  && call_ij)
        {
            // cout << j << endl;
            
            double dXP = wrap_vec(this->X, particles[j].X, Lx);
            double dYP = wrap_vec(this->Y, particles[j].Y, Ly);
            
            double Rij[] = {dXP, dYP};
            double Ui[] = {cos(this->theta), sin(this->theta)};
            double Uj[] = {cos(particles[j].theta), sin(particles[j].theta)};
            double XLiD2 = this->l / 2.0;
            double XLjD2 = particles[j].l / 2.0;
            
            double Rij2 = Rij[0]*Rij[0] + Rij[1]*Rij[1];
            
            double RijEUi = Rij[0] * Ui[0] + Rij[1] * Ui[1];
            double RijEUj = Rij[0] * Uj[0] + Rij[1] * Uj[1];
            double UiEUj = Ui[0] * Uj[0] + Ui[1] * Uj[1];
            double CC = 1 - UiEUj*UiEUj;
            
            double XLanda = 0;
            double XMu = 0;
            double Ro2;
            double dist_ij = 0;
            
            bool brk = false;
            
            while (brk == false)
            {
                // check to see if rods are parallel
                if (CC < exp(-6))
                {
                    if (sqrt(Rij2) >= (XLiD2 + XLjD2)) //rods are parallel and in end-on interaction
                    {
                        XLanda = XLiD2 * sign(RijEUi);
                        XMu = XLjD2 * -sign(RijEUj);
                        
                        if (abs(XMu) > XLjD2){XMu = XLjD2 * sign(XMu);}
                        
                        Ro2 = Rij2 + XLanda*XLanda + XMu*XMu - 2 * XLanda * XMu * UiEUj \
                        + 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
                        dist_ij = sqrt(Ro2);
                        
                        
                        brk = true;
                        break;
                        
                    }
                    //else, rods are parallel and in side-on interaction
                    double XLanda1 = RijEUi + XLjD2 * sign(UiEUj);
                    double XLanda2 = RijEUi - XLjD2 * sign(UiEUj);
                    
                    double XMu1 = -RijEUj + XLiD2 * sign(UiEUj);
                    double XMu2 = -RijEUj - XLiD2 * sign(UiEUj);
                    
                    if (abs(XMu1) > XLjD2) {XMu1 = XLjD2 * sign(XMu1);}
                    if (abs(XMu2) > XLjD2) {XMu2 = XLjD2 * sign(XMu2);}
                    if (abs(XLanda1) > XLiD2) { XLanda1 = XLiD2 * sign(XLanda1);}
                    if (abs(XLanda2) > XLiD2) { XLanda2 = XLiD2 * sign(XLanda2);}
                    
                    XLanda = (XLanda1 + XLanda2) / 2;
                    XMu = (XMu1 + XMu2) / 2;
                    
                    Ro2 = Rij2 + XLanda*XLanda + XMu*XMu - 2 * XLanda * XMu * UiEUj \
                    + 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
                    dist_ij = sqrt(Ro2);
                    brk = true;
                    break;
                }
                
                //cout << "not parallel" << endl;
                //step one after checking if rods are parallel
                XLanda = (RijEUi - UiEUj * RijEUj) / CC;
                XMu = (-RijEUj + UiEUj * RijEUi) / CC;
                
                // are the initially calculated contact points within the bodies of the particles?
                if ((abs(XLanda) <= XLiD2) && (abs(XMu) <= XLjD2))
                {
                    Ro2 = Rij2 + XLanda*XLanda + XMu*XMu - 2 * XLanda * XMu * UiEUj \
                    + 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
                    dist_ij = sqrt(Ro2);
                    brk = true;
                    break;
                }
                
                //cout << "going to aux" << endl;
                
                double AuxIi = abs(XLanda) - XLiD2;
                double AuxIj = abs(XMu) - XLjD2;
                
                if (AuxIi > AuxIj)
                {
                    XLanda = XLiD2 * sign(XLanda);
                    XMu = XLanda * UiEUj - RijEUj;
                    
                    if (abs(XMu) > XLjD2) { XMu = XLjD2 * sign(XMu); }
                    
                }
                
                else
                {
                    XMu = XLjD2 * sign(XMu);
                    XLanda = XMu * UiEUj + RijEUi;
                    
                    if (abs(XLanda) > XLiD2)
                    {XLanda = XLiD2 * sign(XLanda);}
                }
                
                
                Ro2 = Rij2 + XLanda*XLanda + XMu*XMu - 2 * XLanda * XMu * UiEUj \
                + 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
                dist_ij = sqrt(Ro2);
                
                brk = true;
                break;
                
            }
            
            dist_ij = round(dist_ij*10000.0) / 10000.0;
            
            
            if (dist_ij < 2*r0)//eliminates particles from hood that are far away
            {
                new_hood.push_back(particles[j].ID);
            }
            
            //if the particles are in contact, calculate repulsion between them
            if (dist_ij != 0.0 && dist_ij<r0 )
            {
                
                this->calls[j] = false;
                particles[j].calls[i] = false;
                
                double M = -pow(r0, 12.0) / pow(dist_ij, 13.0) + pow(r0, 6.0)/ pow(dist_ij, 7.0);
                
                //contact points
                double XCi = this->X + XLanda * cos(this->theta);
                double YCi = this->Y + XLanda * sin(this->theta);
                double XCj = particles[j].X + XMu * cos(particles[j].theta);
                double YCj = particles[j].Y + XMu * sin(particles[j].theta);
                
                double dXC = wrap_vec(XCi, XCj, Lx);
                double dYC = wrap_vec(YCi, YCj, Ly);
                double UCij[] = {dXC/dist_ij, dYC/dist_ij};
                double U_perp_i_N = Ui[0] * UCij[1] - Ui[1] * UCij[0];
                double U_perp_j_N = Uj[0] * UCij[1] - Uj[1] * UCij[0];
                
                double F_N = 12 * Eo * M;
                
                if (abs(F_N) > F_max){F_N = F_max * sign(F_N);}
                
                this->F_x += F_N * UCij[0];
                this->F_y += F_N * UCij[1];
                particles[j].F_x += -F_N * UCij[0];
                particles[j].F_y += -F_N * UCij[1];
                
                this->torque += F_N * U_perp_i_N * XLanda;
                
                particles[j].torque += -F_N * U_perp_j_N * XMu;
            }
        }
    }
    
    hood = new_hood;
}


void Particle::find_hood_generic(const BINS & b, vector <int> & h, int xi, int yi)

{
    int Xp1 = xi + 1;
    int Xm1 = xi - 1;
    int Yp1 = yi + 1;
    int Ym1 = yi - 1;
    
    if (Xp1 > XI_max){Xp1 = 0;}
    if (Xm1 < 0){Xm1 = XI_max;}
    if (Yp1 > YI_max){Yp1 = 0;}
    if (Ym1 < 0){Ym1 = YI_max;}
    
    h.resize(0);
    
    h.insert(h.end(), b.bins[Ym1][Xm1].begin(), b.bins[Ym1][Xm1].end());
    h.insert(h.end(), b.bins[Ym1][xi].begin(), b.bins[Ym1][xi].end());
    h.insert(h.end(), b.bins[Ym1][Xp1].begin(), b.bins[Ym1][Xp1].end());
    h.insert(h.end(), b.bins[yi][Xm1].begin(), b.bins[yi][Xm1].end());
    h.insert(h.end(), b.bins[yi][xi].begin(), b.bins[yi][xi].end());
    h.insert(h.end(), b.bins[yi][Xp1].begin(), b.bins[yi][Xp1].end());
    h.insert(h.end(), b.bins[Yp1][Xm1].begin(), b.bins[Yp1][Xm1].end());
    h.insert(h.end(), b.bins[Yp1][xi].begin(), b.bins[Yp1][xi].end());
    h.insert(h.end(), b.bins[Yp1][Xp1].begin(), b.bins[Yp1][Xp1].end());
}
