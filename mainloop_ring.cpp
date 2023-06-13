#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 5
#define N2 10000
// #define N 8

double vars[N];
double oldvars[N];
double ktphi[2];
double oldktphi[2];
double chi;
double zeta;
double h;
double b,d;
double errtolmin, errtolmax;
double starting_r, kerr_starting_r;
int zcheck, zinit = 0;

double KT;
int chKT;
int no_isco=0;

double rdisc[N2], emis_del[N2], emis_delnr[N2], inc_del[N2], rmean[N2], rmeannr[N2], intens0[N2], intensnr0[N2];
// double t1,t2,t3;

double eps3, a13, a22, a52, g[4][4], dg_dr[4][4];
int imax = 100;











// #include "eqkerr.h"
// #include "christoffel.h"
#include "eqjohannsen.h"
// #include "eqdCS.h"
// #include "metric.h"
#include "rk45.h"
#include "johannsen.h"
#include "geomet.h"
//#include "gauleg.h"

int main(int argc, char* argv[])
{
    clock_t begin = clock();
    
    
    double chi2, r2;// alpha, beta;
    double cs2, s2;
    double g_tt, g_tp, g_rr, g_thth, g_pp;
    double g_tt_r, g_tp_r, g_pp_r, rhit;
    double E;
    double kt, kphi;
    double scale, oldR;
    double z, oldz, oz2;
    double target_r, hitted_r, hitted_r_old = 1000.0;
    int check, k, target, ec;
    char filename[128];
    int stuck, horizon;
    double deviation;
    FILE *out;
    
    FILE *trajec;
    
    //sprintf(filename, "trajec.txt");
    //trajec = fopen(filename, "w");
    
    double height, delta;
    double I, A, LorentzF, Omega, g_lp;
    // Intensity, proper area, Lorentz factor, orbital velocity of disk, redshift factor, disk's radial thickness
    
    double Emis_delta[2],IIInr[100], ED[100], HD[100];
    
    
    
    double robs[imax], unorm, denom2, johannsen_isco;
    double frac,sign;
    int m, i, j;
    
    M = 1.;
    
    double THETA[imax];
    double RR[imax];
    double frac_glp1, frac_glp2;
    double theta_hitted_array[N2];
    double emissivity[imax];
    double LF[imax];
    double AREA[imax];
    double g_lp_array[imax];
    double III_new[imax];
    double IIInr_new[imax];
    double num_flux[imax];
    
    double III[imax];
    double em_profile[imax];
    double em_profile_new[imax];
    double num_flux_new[imax];
    double int_prof[imax];
    double int_prof_new[imax];
    double robs_gauleg[imax];
    double r_ring;
    
    
    double hx_start, hy_start, hy_end, hx_end, Delhy;
    double rr_start, theta_start;
    double omega_cr;
    double omega_glp;
    double omega_cr_term;
    double beta;
    
    double e_phit;
    double denom_e_phit;
    double denom_fact1_e_phit;
    double denom_fact2_e_phit;
    double denom_fact3_e_phit;
    double numeretor_e_phit;
    double sigma;
    double omega_small;
    
    double p1;
    double nume_fact1_ephi1;
    double nume_fact2_ephi1;
    double nume_fact3_ephi1;
    double nume_fact4_ephi1;
    double nume_ephi1;
    double denom_fact1_e_phi1;
    double denom_fact2_e_phi1;
    double denom_fact3_e_phi1;
    double denom_e_phi1;
    double e_phi1;
    
    double rho_sq; // this is the same factor as in johansen meteric Sigma = r*r + a*a*cos(theta)*cos(theta)
    double f;
    double A1;
    double A2;
    double A5;
    double B;
    double rhotilda_sq;
    double a_sq;
    double e2psi ;
    double e2nu ;
    
    
    
    for(int j = 0; j<imax; j = j+1){
        III[j] = 0.0;
        III_new[j] = 0.0;
        num_flux[j] = 0.0;
        //num_flux_new[j] = 0.0;
        AREA[j] = 0.0;
        LF[j] = 0.0;
        em_profile[j] = 0.0;
        //em_profile_new[j] = 0.0;
        int_prof[j] = 0.0;
        //int_prof_new[j] = 0.0;
    }
    
    
    
    
    
    /*-------------computing (r, theta) for coronal geometery--------------- */
    
    
    
    
    a = atof(argv[1]);      /* spin parameter */
    
    
    chi = a;
    
    a_sq = a*a;
    
    // /* Deformation Parameters */
    
    a13 =atof(argv[2]);
    
    a22 = 0.0;//atof(argv[5]);
    a52 = 0.0;//atof(argv[6]);
    eps3 = 0.0;//atof(argv[7]);
    
    
    
    find_isco(50., johannsen_isco);
    while (no_isco) {
        no_isco = 0;
        if (a13+a22+a52+eps3 < 0) sign = -1.0;
        else if (a13+a22+a52+eps3 > 0) sign = 1.0;
        else sign = 0.0;
        a13 -= sign*a13/1000.0;
        a22 -= sign*a22/1000.0;
        a52 -= sign*a52/1000.0;
        eps3 -= sign*eps3/1000.0;
        find_isco(50., johannsen_isco);
        
    }
    
    kerr_starting_r = johannsen_isco;
    
    
    
    //fprintf(out, "%f\n", johannsen_isco);
    printf(" a = %f, a13=%f, a22=%f, a52=%f, eps3=%f\n", a, a13, a22, a52, eps3);
    
    starting_r = johannsen_isco * kerr_starting_r / kerr_rms(a);
    
    printf(" starting radii=%f, isco=%f, horizon=%f\n", starting_r, johannsen_isco, 1.0 + sqrt(1.0 - a*a));
    
    
    
    //gauleg(kerr_starting_r, 1000, robs_gauleg);
    
    
    for (i = 0; i <= imax - 1; i++) {
        robs[i] = pow((double)i/99.,3) * (1000. - kerr_starting_r) + kerr_starting_r;
    }
    
    
    /*----------computing Area and lorentz factor arrays---------------*/
    
    for(i = 0; i<imax; i++){
        
        metric(robs[i], M_PI/2.0, g);
        metric_rderivatives(robs[i], M_PI/2.0, dg_dr);
        
        
        AREA[i] = 2 * M_PI * sqrt(g[1][1] * g[3][3])*(robs[i+1] - robs[i]);
        
        Omega = (-dg_dr[0][3] + sqrt(dg_dr[0][3]*dg_dr[0][3] - dg_dr[0][0]*dg_dr[3][3])) / dg_dr[3][3];
        
        LF[i] = pow(pow(Omega + g[0][3] / g[3][3], 2) * g[3][3] * g[3][3] / (g[0][0]*g[3][3] - g[0][3]*g[0][3]) + 1.0, -0.5);
        
    }
    
    /*-----------------------------------------------------------------------------*/
    
    double rh;
    rh = 1.0 + sqrt(1.0 - a*a);
    int h_index;
    
    for(h_index = 0; h_index<=33; h_index++){
        
        hy_start = pow(h_index/249.,2) * (500. - rh*1.75 ) + rh*1.75;

        
        for(r_ring = 1; r_ring<=1; r_ring += 25./15.){
            
            for(int j = 0; j<imax; j = j+1){
                num_flux[j] = 0.0;
                III[j] = 0.0;

                em_profile[j] = 0.0;

                int_prof[j] = 0.0;

            }
            
            sprintf(filename,"a%.05e.h_%.03e_r_%.03e_e_%.03e.a13_%.03e.a22_%.03e.a52_%.03e.dat",a,hy_start,r_ring,eps3,a13,a22,a52);
            
            out = fopen(filename,"w");
            
            for(hx_start = r_ring - 0.25; hx_start<=r_ring + 0.25; hx_start += 0.25){ // radial distribution of emitting points on the coronal disk
                
                
                geomet(hx_start, hy_start, rr_start, theta_start);
                
                printf("hx_start %f, r_ring %f ,hy_start %f, rr_start %f, theta_start %f\n", hx_start, r_ring ,hy_start, rr_start, theta_start);
                
                
                
                for(beta = 0.0000001; beta<6.2831853; beta += 6.2831853/15.){ // the loop on the beta
                    
                    double cos_alpha;
                    
                    for(cos_alpha = -1; cos_alpha< 1.; cos_alpha += 1./25999.){ // loop on delta
                        
                        delta = acos(cos_alpha);
                        
                        
                        chi2 = chi*chi;
                        
                        h = -1.0;
                        errtolmin = 1.0e-6;
                        errtolmax = 1.0e-8;
                        
                        d = 2000.0;
                        
                        scale = 1.;
                        check = 0;
                        target = 0;
                        
                        //  Initial values
                        
                        horizon = 0;
                        zinit = 0;  //  additional parameter for runge-kutta
                        ec = 0;     //
                        
                        vars[0] = rr_start; //r-cordinate
                        if(theta_start == 0.0){vars[1] = 1E-8;}
                        else{vars[1] = theta_start;}//asin(0.7747764 / height); theta-cordinate}
                        vars[4] = 0.0; //phi-cordinate
                        
                        r2 = vars[0]*vars[0];
                        
                        /*---computing some terms in Johansen metric------*/
                        
                        denom2 = sqrt(r2 - 2*vars[0] + a*a); // this term is basically sqrt(DELTA)
                        rho_sq = r2 + a_sq*cos(vars[1])*cos(vars[1]) ;
                        
                        f = eps3/vars[0];
                        A1 = 1 + (a13/vars[0]);
                        A2 = 1 + (a22/vars[0]);
                        A5 = 1 + (a52/vars[0]);
                        B = ((r2 + a_sq)*A1) - (a_sq*A2*sin(vars[1])*sin(vars[1]));
                        rhotilda_sq = rho_sq + f;
                        /*---------------------------------------------*/
                        
                        /*---r and theta velociites of photons--------*/
                        vars[2] = (cos(delta)*denom2) / sqrt(rhotilda_sq); //r velocity
                        vars[3] = sin(delta)*cos(beta)/ sqrt(rhotilda_sq); //theta velocity
                        /*-----------------------------------------------*/
                        
                        
                        /*---computing kphi for non-kerr case ------*/
                        
                        omega_small = (a*( ( (r2 + a_sq)*A1*A2) - (denom2*denom2) ) ) / ( ((r2 + a_sq)*(r2 + a_sq)*A1*A1) - (a_sq*denom2*denom2*sin(vars[1])*sin(vars[1])));
                        
                        e2psi = ((((r2 + a_sq)*(r2 + a_sq)*A1*A1) - (a_sq*denom2*denom2*sin(vars[1])*sin(vars[1])))*rhotilda_sq*sin(vars[1])*sin(vars[1])) / (B*B);
                        
                        e2nu = (e2psi*omega_small*omega_small) + ( rhotilda_sq*( (denom2*denom2) - (a_sq*A2*A2*sin(vars[1])*sin(vars[1])) ) / (B*B) );
                        
                        
                        
                        
                        omega_cr_term = pow(hx_start, 1.5);
                        
                        omega_cr =0.0;// 1.0/(fabs(a) + omega_cr_term);
                        
                        
                        numeretor_e_phit = omega_cr*(1/sqrt(e2nu)) ; // e^-nu = 1/sqrt(e^2nu)
                        
                        
                        denom_fact1_e_phit = (omega_cr - omega_small )*(omega_cr - omega_small );
                        
                        denom_fact2_e_phit = e2psi*(1/e2nu);
                        
                    
                        denom_e_phit = sqrt(1 - (denom_fact2_e_phit*denom_fact1_e_phit));
                        
                        
                        e_phit = numeretor_e_phit / denom_e_phit ;
                        
                        p1 = sin(delta)*sin(beta);
                        
                        nume_fact1_ephi1 = (1/sqrt(fabs(e2nu)))*(1/sqrt(fabs(e2psi)));
                        nume_fact2_ephi1 = e2nu;
                        nume_fact3_ephi1 = omega_cr*omega_small*e2psi;
                        nume_fact4_ephi1 = omega_small*omega_small*e2psi;
                        
                        nume_ephi1 = nume_fact1_ephi1*(nume_fact2_ephi1 + nume_fact3_ephi1 - nume_fact4_ephi1 );
                        
                        denom_fact1_e_phi1 = e2nu ;
                        denom_fact2_e_phi1 = e2psi;
                        denom_fact3_e_phi1 = (omega_cr - omega_small)*(omega_cr - omega_small) ;
                        
                        denom_e_phi1 = sqrt(denom_fact1_e_phi1 - denom_fact2_e_phi1*denom_fact3_e_phi1);
                        
                        e_phi1 = nume_ephi1 / denom_e_phi1 ;
                        
                        kphi = e_phit + p1*e_phi1;
                        
                        /*-------------------------------------------------------------------------------------*/
                        
                        z = vars[0]*cos(vars[1]);
                        oldz = 10000.;
                        oz2 = oldz;
                        
                        cs2 = cos(vars[1])*cos(vars[1]);
                        s2 = sin(vars[1])*sin(vars[1]);
                        
                        metric(vars[0], vars[1], g);
                        metric_rderivatives(vars[0], vars[1], dg_dr);
                        
                        g_thth = g[2][2];
                        g_tt = g[0][0];
                        g_tp = g[0][3];
                        g_rr = g[1][1];
                        g_pp = g[3][3];
                        
                        E = sqrt(g_tp*g_tp*kphi*kphi - g_tt*(g_rr*vars[2]*vars[2]+g_thth*vars[3]*vars[3]+g_pp*kphi*kphi));
                        
                        kt = -(g_tp*kphi+E)/g_tt;
                        
                        b = -(g_pp*kphi+g_tp*kt)/(g_tt*kt+g_tp*kphi);
                        
                        ktphi[0] = kt/E;
                        ktphi[1] = kphi/E;
                        
                        
                        unorm = -(g_rr * pow(vars[2], 2) + g_thth * pow(vars[3], 2) + g_pp * pow(kphi, 2) + 2.0 * g_tp * kt * kphi) / (g_tt * pow(kt, 2));
                        
                    
                        
                        vars[2] /= E;
                        vars[3] /= E;
                        

                        
                        zcheck = 0;   // additional parameter for runge-kutta
                        
                                                
                        for(k = 0; !horizon && z*oldz >= 0. && vars[0] <= d && vars[0] > 1.0 + sqrt(1.0 - a*a) + 0.001; k++) {
                            
                            cache();
                            rk();
                            
                      
                            oldz = z;
                            
                            z = vars[0]*cos(vars[1]);
                            if (z == oldz) {
                                if (z == oz2) {
                                    ec++;  // for case, when photon goes into horizon
                                }
                                else
                                {
                                    oz2 = z;
                                    ec = 0;
                                }
                            }
                            
                            
                            if (ec > 5*pow(10, 3)) {    //  horizon case, !!! 9 th power is the temporary solution
                                printf(" !!! horizon = %f\n", vars[0]);
                                horizon = 1;
                                
                            }
                            
                            if (z*oldz <= 0. && fabs(z - oldz) > 1E-8) {
                                zcheck = 1;
                                zinit = 1;
                                rcache();
                                z = oldz;
                            }
                        }
                        
                        //abort();
                        
                        if (z * oldz <= 0. && vars[0] >= kerr_starting_r*0.9 && vars[0] <= 1000.) {
                            
                            rhit = find_radii(vars, oldvars);
                            zcheck = 0;

                            
                            //-----------computing g_lp-------------------
                            
                            
                            
                            
                            metric(rr_start, theta_start, g);
                            
                            //omega_cr_term = pow(hx_start, 1.5);
                            
                            //omega_cr = 1.0/(a + omega_cr_term);
                            
                            frac_glp1 = g[0][0] + (2*g[0][3]*omega_cr) + (g[3][3]*omega_cr*omega_cr);
                            
                            
                            metric_rderivatives(rhit, M_PI/2.0, dg_dr);
                            Omega = (-dg_dr[0][3] + sqrt(dg_dr[0][3]*dg_dr[0][3] - dg_dr[0][0]*dg_dr[3][3])) / dg_dr[3][3];
                            metric(rhit, M_PI/2., g);
                            frac_glp2 = g[0][0] + (2*g[0][3]*Omega) + (g[3][3]*Omega*Omega);
                            g_lp = sqrt(frac_glp1 / frac_glp2);
                            

                            
                            for(i = 0; i<=imax - 2; i++){
                                
                                if(robs[i] < rhit && robs[i+1] > rhit)
                                {
                                    num_flux[i] = num_flux[i] + 1.*g_lp*g_lp;
                                    III[i] = III[i] + 1.;
                                }
                            }
                        }
                    }
                }
            }
            
            /*----interpolation-------*/

            
           
            
            for(i = 0; i<imax; i++){
                
                em_profile[i] = (num_flux[i])/(AREA[i]*LF[i]);
                
                int_prof[i] = (III[i])/(AREA[i]*LF[i]);
                
            }
            
            for(i = 0; i< imax; i++){
                
                
                fprintf(out,"%f %f %f\n",robs[i], em_profile[i], int_prof[i]);
                
            }
            
            clock_t end = clock();
            double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            printf(" Execution time is %f s\n", time_spent);
            
            fprintf(out, "\n");
            fclose(out);
            
        }  // closing the loop over hx_end (r) i.e the radius of the ring
        
    } // closing the loop over hy_start (h) i.e height of the ring
    return 0;
    
}


//////////////////////////////////////////////////
/////////////////  END OF MAIN ///////////////////
/////////////////////////////////////////////////
