// #include "metric.h"

const  double Pi  = 3.141592653589793;
const  double Piby2 = 1.5707963267948966192;
double a, M;



void cache();
void rcache();

double R(double vars[5]);

double kerr_rms(double a);
double find_hitting_angle(double nvars[5], double ovars[5]);
double find_radii(double nvars[5], double ovars[5]);
double interp_lin_1d(double ifac_r, double rlo, double rhi);
double cal_q(double delta, double h);
double cal_delta_inc(double delta, double theta, double h, double r);
double cal_kt(double r, double theta);
double cal_gi(double delta, double theta, double h, double r);
double r0_grid(int i);
double r_mean(int i);
double I0(int i);
double II(double r);
double flat_I(double r,  double h);
void sortd(double ds[]);



// extern double a13;

// void metric( double spin,  double spin2,  double eps3,  double a13,  double a22,  double a52,  double z1,  double z2,  double mn[][4]);
// void metric_rderivatives( double spin,  double spin2,  double eps3,  double a13,  double a22,  double a52,  double z1,  double z2,  double dmn[][4]);
// void metric(double z1, double z2, double mn[][4]);
// void metric_rderivatives(double z1, double z2, double dmn[][4]);
void metric(double r, double th, double mn[][4]);
void metric_rderivatives(double r, double th, double dmn[][4]);
void find_isco(double z1, double& isco);

// void find_isco( double spin,  double spin2,  double eps3,  double a13,  double a22,  double a52,  double z1,  double& isco);



double gtt(double a, double r, double x, double M);
double grr(double a, double r, double x, double M);
double gxx(double a, double r, double x, double M);
double gpp(double a, double r, double x, double M);
double gtp(double a, double r, double x, double M);

double gltt(double a, double r, double x, double M);
double glrr(double a, double r, double x, double M);
double glxx(double a, double r, double x, double M);
double glpp(double a, double r, double x, double M);
double gltp(double a, double r, double x, double M);

double dgltt_dr(double a, double r, double x, double M);
double dglrr_dr(double a, double r, double x, double M);
double dglxx_dr(double a, double r, double x, double M);
double dglpp_dr(double a, double r, double x, double M);
double dgltp_dr(double a, double r, double x, double M);

double dgltt_dx(double a, double r, double x, double M);
double dglrr_dx(double a, double r, double x, double M);
double dglxx_dx(double a, double r, double x, double M);
double dglpp_dx(double a, double r, double x, double M);
double dgltp_dx(double a, double r, double x, double M);




double gtt(double a, double r, double x, double M) {
  return(-(pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2)) + 
     pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2))/
   (pow(r,5)*(pow(a,2) + (-2 + r)*r)*
     (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2)));
}

double grr(double a, double r, double x, double M) {
  return ((pow(a,2) + (-2 + r)*r)*(a52 + pow(r,2)))/
   (r*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2)));
}

double gxx(double a, double r, double x, double M) {
  return r/(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2));
}

double gpp(double a, double r, double x, double M) {
  return (pow(r,4)*(pow(a,2) + (-2 + r)*r) - 
     pow(a,2)*pow(a22 + pow(r,2),2)*pow(sin(x),2))/
   (pow(r,3)*(pow(a,2) + (-2 + r)*r)*
     (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2));
}

double gtp(double a, double r, double x, double M) {
  return -((a*(pow(r,5)*(a22 + 2*r) + a13*pow(r,2)*(a22 + pow(r,2)) + 
         pow(a,2)*(a22*pow(r,3) + a13*(a22 + pow(r,2)))))/
     (pow(r,4)*(pow(a,2) + (-2 + r)*r)*
       (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))));
}




double gltt(double a, double r, double x, double M) {
  return -((r*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*
       (pow(r,4)*(pow(a,2) + (-2 + r)*r) - 
         pow(a,2)*pow(a22 + pow(r,2),2)*pow(sin(x),2)))/
     pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
       pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2));
}

double glrr(double a, double r, double x, double M) {
  return (r*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2)))/
   ((pow(a,2) + (-2 + r)*r)*(a52 + pow(r,2)));
}

double glxx(double a, double r, double x, double M) {
  return eps3/r + pow(r,2) + pow(a,2)*pow(cos(x),2);
}

double glpp(double a, double r, double x, double M) {
  return ((eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2)*
     (pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2) - 
       pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2)))/
   (r*pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
       pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2));
}

double gltp(double a, double r, double x, double M) {
  return -((a*(pow(r,5)*(a22 + 2*r) + a13*pow(r,2)*(a22 + pow(r,2)) + 
         pow(a,2)*(a22*pow(r,3) + a13*(a22 + pow(r,2))))*
       (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2))/
     pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
       pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2));
}




double dgltt_dr(double a, double r, double x, double M) {
  return -((r*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*
        (pow(r,4)*(-2 + 2*r) + 4*pow(r,3)*(pow(a,2) + (-2 + r)*r) - 
          4*pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2)))/
      pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
        pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2)) + 
   (2*r*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*
      (3*pow(r,2)*(pow(a,2) + pow(r,2)) + 2*r*(a13 + pow(r,3)) - 
        2*pow(a,2)*pow(r,2)*pow(sin(x),2) - 
        pow(a,2)*(a22 + pow(r,2))*pow(sin(x),2))*
      (pow(r,4)*(pow(a,2) + (-2 + r)*r) - 
        pow(a,2)*pow(a22 + pow(r,2),2)*pow(sin(x),2)))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),3) - 
   (r*(3*pow(r,2) + pow(a,2)*pow(cos(x),2))*
      (pow(r,4)*(pow(a,2) + (-2 + r)*r) - 
        pow(a,2)*pow(a22 + pow(r,2),2)*pow(sin(x),2)))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2) - 
   ((eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*
      (pow(r,4)*(pow(a,2) + (-2 + r)*r) - 
        pow(a,2)*pow(a22 + pow(r,2),2)*pow(sin(x),2)))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2);
}

double dglrr_dr(double a, double r, double x, double M) {
  return (r*(3*pow(r,2) + pow(a,2)*pow(cos(x),2)))/
    ((pow(a,2) + (-2 + r)*r)*(a52 + pow(r,2))) - 
   (2*pow(r,2)*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2)))/
    ((pow(a,2) + (-2 + r)*r)*pow(a52 + pow(r,2),2)) - 
   (r*(-2 + 2*r)*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2)))/
    (pow(pow(a,2) + (-2 + r)*r,2)*(a52 + pow(r,2))) + 
   (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))/
    ((pow(a,2) + (-2 + r)*r)*(a52 + pow(r,2)));
}

double dglxx_dr(double a, double r, double x, double M) {
  return -(eps3/pow(r,2)) + 2*r;
}

double dglpp_dr(double a, double r, double x, double M) {
  return (-2*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2)*
      (pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2) - 
        pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2))*
      (3*pow(r,2)*(pow(a,2) + pow(r,2)) + 2*r*(a13 + pow(r,3)) - 
        2*pow(a,2)*pow(r,2)*pow(sin(x),2) - 
        pow(a,2)*(a22 + pow(r,2))*pow(sin(x),2)))/
    (r*pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
        pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),3)) + 
   ((eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2)*
      (6*pow(r,2)*pow(pow(a,2) + pow(r,2),2)*(a13 + pow(r,3)) + 
        4*r*(pow(a,2) + pow(r,2))*pow(a13 + pow(r,3),2) - 
        pow(a,2)*pow(r,6)*(-2 + 2*r)*pow(sin(x),2) - 
        6*pow(a,2)*pow(r,5)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2)))/
    (r*pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
        pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2)) + 
   ((3*pow(r,2) + pow(a,2)*pow(cos(x),2))*pow(sin(x),2)*
      (pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2) - 
        pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2)))/
    (r*pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
        pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2)) - 
   ((eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2)*
      (pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2) - 
        pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2)))/
    (pow(r,2)*pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
        pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2));
}

double dgltp_dr(double a, double r, double x, double M) {
  return (2*a*(pow(r,5)*(a22 + 2*r) + a13*pow(r,2)*(a22 + pow(r,2)) + 
        pow(a,2)*(a22*pow(r,3) + a13*(a22 + pow(r,2))))*
      (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2)*
      (3*pow(r,2)*(pow(a,2) + pow(r,2)) + 2*r*(a13 + pow(r,3)) - 
        2*pow(a,2)*pow(r,2)*pow(sin(x),2) - 
        pow(a,2)*(a22 + pow(r,2))*pow(sin(x),2)))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),3) - 
   (a*(pow(r,5)*(a22 + 2*r) + a13*pow(r,2)*(a22 + pow(r,2)) + 
        pow(a,2)*(a22*pow(r,3) + a13*(a22 + pow(r,2))))*
      (3*pow(r,2) + pow(a,2)*pow(cos(x),2))*pow(sin(x),2))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2) - 
   (a*(2*a13*pow(r,3) + 2*pow(r,5) + 5*pow(r,4)*(a22 + 2*r) + 
        2*a13*r*(a22 + pow(r,2)) + pow(a,2)*(2*a13*r + 3*a22*pow(r,2)))*
      (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),2))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2);
}





double dgltt_dx(double a, double r, double x, double M) {
  return (2*pow(a,2)*pow(r,2)*cos(x)*(pow(r,4)*(pow(a,2) + (-2 + r)*r) - 
        pow(a,2)*pow(a22 + pow(r,2),2)*pow(sin(x),2))*sin(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2) + 
   (2*pow(a,2)*r*pow(a22 + pow(r,2),2)*
      (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*sin(x)*cos(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2) - 
   (4*pow(a,2)*pow(r,2)*(a22 + pow(r,2))*
      (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*sin(x)*
      (pow(r,4)*(pow(a,2) + (-2 + r)*r) - 
        pow(a,2)*pow(a22 + pow(r,2),2)*pow(sin(x),2))*cos(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),3);
}

double dglrr_dx(double a, double r, double x, double M) {
  return (-2*pow(a,2)*pow(r,2)*cos(x)*sin(x))/((pow(a,2) + (-2 + r)*r)*(a52 + pow(r,2)));
}

double dglxx_dx(double a, double r, double x, double M) {
  return -2*pow(a,2)*cos(x)*sin(x);
}

double dglpp_dx(double a, double r, double x, double M) {
  return (-2*pow(a,2)*cos(x)*pow(sin(x),2)*
      (pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2) - 
        pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2))*sin(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2) + 
   (4*pow(a,2)*(a22 + pow(r,2))*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*
      pow(sin(x),3)*(pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2) - 
        pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2))*
      cos(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),3) - 
   (2*pow(a,2)*pow(r,5)*(pow(a,2) + (-2 + r)*r)*
      (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),3)*
      cos(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2) + 
   (2*(eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*sin(x)*
      (pow(pow(a,2) + pow(r,2),2)*pow(a13 + pow(r,3),2) - 
        pow(a,2)*pow(r,6)*(pow(a,2) + (-2 + r)*r)*pow(sin(x),2))*
      cos(x))/
    (r*pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
        pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2));
}

double dgltp_dx(double a, double r, double x, double M) {
  return (2*pow(a,3)*r*(pow(r,5)*(a22 + 2*r) + a13*pow(r,2)*(a22 + pow(r,2)) + 
        pow(a,2)*(a22*pow(r,3) + a13*(a22 + pow(r,2))))*cos(x)*pow(sin(x),2)*
      sin(x))/pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2) - 
   (4*pow(a,3)*r*(a22 + pow(r,2))*
      (pow(r,5)*(a22 + 2*r) + a13*pow(r,2)*(a22 + pow(r,2)) + 
        pow(a,2)*(a22*pow(r,3) + a13*(a22 + pow(r,2))))*
      (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*pow(sin(x),3)*
      cos(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),3) - 
   (2*a*(pow(r,5)*(a22 + 2*r) + a13*pow(r,2)*(a22 + pow(r,2)) + 
        pow(a,2)*(a22*pow(r,3) + a13*(a22 + pow(r,2))))*
      (eps3 + pow(r,3) + pow(a,2)*r*pow(cos(x),2))*sin(x)*cos(x))/
    pow((pow(a,2) + pow(r,2))*(a13 + pow(r,3)) - 
      pow(a,2)*r*(a22 + pow(r,2))*pow(sin(x),2),2);
}


void metric(double r, double th, double mn[][4])
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
  double g_tt, g_rr, g_thth, g_pp, g_tp;
  double spin = chi;
  
  t1 = cos(th);
  t2 = spin * spin;
  t3 = r * r;
  t4 = pow(t3, 0.2e1);
  t5 = t3 * t4;
  t1 = t2 * pow(t1, 0.2e1);
  t6 = (t1 + t3) * r + eps3;
  t7 = a22 + t3;
  t8 = sin(th);
  t9 = t3 + t2;
  t8 = pow(t8, 0.2e1);
  t10 = r * t3 + a13;
  t11 = -t2 * r * t7 * t8 + t10 * t9;
  t12 = -0.2e1 * r + t9;
  t13 = a52 + t3;
  t14 = 0.1e1 / r;
  t15 = t2 * a22;
  t16 = 0.1e1 / t12;
  t13 = 0.1e1 / t13;
  t11 = pow(t11, -0.2e1);
  
  g_tt = t6 * (t2 * pow(t7, 0.2e1) * t8 + 0.2e1 * r * t4 - t4 * t9) * r * t11;
  g_rr = t6 * r * t16 * t13;
  g_thth = t14 * eps3 + t1 + t3;
  g_pp = (pow(t9, 0.2e1) * pow(t10, 0.2e1) - t2 * t5 * t12 * t8) * t6 * t8 * t11 * t14;
  g_tp = -(0.2e1 * t5 + t3 * (a13 * (a22 + t2) + ((r * a22 + a13) * r + t15) * r) + t15 * a13) * spin * t6 * t8 * t11;

  mn[0][0] = g_tt;
  mn[0][3] = g_tp;
  mn[1][1] = g_rr;
  mn[2][2] = g_thth;
  mn[3][0] = mn[0][3];
  mn[3][3] = g_pp;
}

void metric_rderivatives(double r, double th, double dmn[][4])
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23;
  double dgttdr, dgtpdr, dgppdr;
  double spin = chi;
  
  t1 = cos(th);
  t2 = spin * spin;
  t3 = r * r;
  t4 = pow(t3, 0.2e1);
  t5 = t3 * t4;
  t1 = t2 * pow(t1, 0.2e1);
  t6 = (t1 + t3) * r + eps3;
  t7 = a22 + t3;
  t8 = sin(th);
  t9 = t3 + t2;
  t10 = -0.2e1 * r + t9;
  t8 = pow(t8, 0.2e1);
  t11 = r * t3 + a13;
  t12 = t9 * t11;
  t13 = -t2 * r * t7 * t8 + t12;
  t1 = t1 / 0.3e1 + t3;
  t14 = 0.2e1 / 0.3e1 * t2;
  t15 = 0.3e1 / 0.5e1 * t2;
  t16 = (t15 + t3) * r + 0.2e1 / 0.5e1 * a13;
  t15 = r * t16 - t15 * (a22 / 0.3e1 + t3) * t8;
  t13 = 0.1e1 / t13;
  t17 = pow(t13, 0.2e1);
  t18 = 0.3e1 * t1;
  t19 = a22 + t2;
  t20 = r * a22;
  t21 = t2 * a22;
  t22 = 0.1e1 / 0.2e1;
  t23 = t17 * t8;
  t9 = -t2 * t5 * t10 * t8 + pow(t9, 0.2e1) * pow(t11, 0.2e1);
  t11 = 0.1e1 / r;
  
  dgttdr = t17 * ((t6 * (0.10e2 * r * t15 * t13 - 0.1e1) - t18 * r) * (-t2 * pow(t7, 0.2e1) * t8 + t10 * t4) - 0.6e1 * t6 * (t3 * ((-0.5e1 / 0.3e1 + r) * r + t14) - t14 * t7 * t8) * t3);
  dgtpdr = t23 * spin * ((0.20e2 * t6 * t15 * t13 - 0.6e1 * t1) * (t22 * (t3 * (a13 * t19 + ((a13 + t20) * r + t21) * r) + t21 * a13) + t5) - 0.12e2 * t6 * ((t19 / 0.6e1 + t3 / 0.3e1) * a13 + t4 + t20 * (0.5e1 / 0.12e2 * t3 + t2 / 0.4e1)) * r);
  dgppdr = 0.10e2 * t23 * t6 * (-t15 * t9 * t11 * t13 + t12 * t16 - 0.4e1 / 0.5e1 * t4 * ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t2) * t2 * t8) + t23 * t9 * t11 * (-t11 * t6 + t18);
  
  dmn[0][0] = dgttdr;
  dmn[0][3] = dgtpdr;
  dmn[3][0] = dmn[0][3];
  dmn[3][3] = dgppdr;
}

void find_isco(double z1, double& isco) 
{
  // printf(" %f, %f\n", chi, a13);
  
  int i, j, casenum=1, stop=0, count=0;
  double spin = chi;
  double detol = 1.0e-8;
  double rll,rul,rinit=z1,rnew,rold,rstep=1.e-5;
  double mn[4][4],dmn[4][4];
  double mnp[4][4],mnm[4][4],dmnp[4][4],dmnm[4][4];
  double ep,em,eold,enew,omp,omm,omold,omnew;
  double epsq,emsq;
  double deold,denew;
  double dr=1.0e-5;
  double sqspin=spin*spin;

  if(spin>0.)
    rll = 1.+sqrt(1.-sqspin);
  else if(spin<0.)
    rll = 1.+sqrt(1.-sqspin);
  else
    rll = 1.;

  if(a13>-5.7){
    rul = rinit; rold = rul;
    metric(rold+dr,Pi/2.,mnp);
    metric(rold-dr,Pi/2.,mnm);
    metric_rderivatives(rold+dr,Pi/2.,dmnp);
    metric_rderivatives(rold-dr,Pi/2.,dmnm);
    omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3]; 
    omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3]; 
    ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(-mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp);
    em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(-mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm);
    deold = 0.5*(ep-em)/dr;
    
    do{
      count++;
      if(count>40){
        printf("No convergence after %i iterations.\n",count);
        no_isco=1;
        break;
      }

      rnew = (rll+rul)/2.;
      metric(rnew+dr,Pi/2.,mnp);
      metric(rnew-dr,Pi/2.,mnm);
      metric_rderivatives(rnew+dr,Pi/2.,dmnp);
      metric_rderivatives(rnew-dr,Pi/2.,dmnm);
      omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3]; 
      omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3]; 
      ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(-mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp);
      em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(-mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm);
      denew = 0.5*(ep-em)/dr;
      
      if(fabs(denew)<fabs(detol)) {
          //printf("denew = %e, deold = %e\n",denew,deold);
          stop = 1;
      }
      else if((denew*deold)>0.0) {
        if(rnew<rold)
          rul = rnew;
        else if(rnew>rold) 
          rll = rnew;
        else
          printf("rold=rnew? rold = %e, rnew = %e",rold,rnew);
      }
      else if((denew*deold)<0.0) {
        if(rnew<rold)
          rll = rnew;
        else if(rnew>rold)
          rul = rnew;
        else
          printf("rold=rnew? rold = %e, rnew = %e",rold,rnew);
      }
      else
        printf("Compare enew and eold. eold = %e, enew = %e\n",deold,denew);
     
       
      rold = rnew;
      omold = omnew;
      deold = denew;

    }while(stop==0);
  }
  else {
    rold = 4.0; // ISCO radius for these special cases is never above this, so we start here. -SN
    
    do{
      rnew = rold-rstep;
      
      if(rnew<rll){
        printf("Couldn't find ISCO? rnew = %e, rll = %e\n",rnew,rll);
        break;
      }
    
      metric(rnew+dr,Pi/2.,mnp);
      metric(rnew-dr,Pi/2.,mnm);
      metric_rderivatives(rnew+dr,Pi/2.,dmnp);
      metric_rderivatives(rnew-dr,Pi/2.,dmnm);
      omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3]; 
      omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3]; 
      epsq = -mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp;
      emsq = -mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm;
      
      if(epsq>0. && emsq>0.){
        ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(epsq);
        em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(emsq);
        denew = 0.5*(ep-em)/dr;
        
        if(fabs(denew)<fabs(detol)) {
          //printf("denew = %e, deold = %e\n",denew,deold);
          stop = 1;
          //break;
        }
        else {
          //printf("%e\n",rnew);
          rold = rnew;
        }
      }
      else{
        printf("epsq = %e, emsq = %e, rnew = %e\n",epsq,emsq,rnew);
        stop = 1;
      }
      
    }while(stop==0);
  }

    
  if(stop==1)
    isco=rnew;
}













double kerr_rms(double a){
  //   accounts for negative spin
  double sign = 1.0;
  if (a<0) {
     sign = -1.0;
  }

  double Z1 = 1.0+pow(1.0-a*a,1.0/3.0)*(pow(1.0+a,1.0/3.0)+pow(1.0-a,1.0/3.0));
  double Z2=sqrt((3.0*a*a)+(Z1*Z1));

  return 3.0+Z2-sign*sqrt((3.0-Z1)*(3.0+Z1+(2*Z2)));
}



double flat_I(double r,  double h) {
  return 1. / pow(pow(r/h, 2.) + 1., 3./2.) / (2. * M_PI * pow(h, 2.));
}


void sortd(double ds[]){
  int i, j;
  double tmp;
  for (i=0;i<3;i++) 
    for (j=i+1; j<3; j++)
      if (ds[i] > ds[j]) {
        tmp = ds[i];
        ds[i] = ds[j];
        ds[j] = tmp;
      }
  
}


double cal_q(double delta, double height) {
  double a = chi;
  // double q = sqrt(a*a - 2.0*h + h*h) * sin(delta) / sqrt(a*a + h*h);
  double q = (a*a + height*height) * sin(delta);
  q *= q;
  q /= (height*height - 2*height + a*a);
  q -= a*a;
  q = sqrt(fabs(q));
  // q = sqrt((h*h - 2*h + a*a) / (a*a + h*h)) * sin(delta);
  return q;
}

double cal_delta_inc(double delta, double theta, double height, double r) {
  // double q = cal_q(delta, height);
  // double kt = cal_kt(r, theta);

  // double E, Lz;
  // double kt1,kp1,kth1,kr1;
  // double th1 = fabs(oldvars[1]), th2 = fabs(vars[1]), coef1, coef2;
  // th1 = M_PI/2.0 - th1;
  // th2 = th2 - M_PI/2.0;

  // coef1 = th1 / (th1+th2);
  // coef2 = th2 / (th1+th2);

  // kt1 = coef1*oldktphi[0] + coef2*ktphi[0];
  // kp1 = coef1*oldktphi[1] + coef2*ktphi[1];
  // kr1 = coef1*oldvars[2] + coef2*vars[2];
  // kth1 = coef1*oldvars[3] + coef2*vars[3];

  // kt1 = oldktphi[0];
  // kp1 = oldktphi[1];
  // kr1 = oldvars[2];
  // kth1 = oldvars[3];


  metric(r, theta, g);
  metric_rderivatives(r, theta, dg_dr);
  // double E1 = sqrt(g[0][3]*g[0][3]*kp1*kp1 - g[0][0]*(g[1][1]*kr1*kr1+g[2][2]*kth1*kth1+g[3][3]*kp1*kp1));

  double Omega = (-dg_dr[0][3] + sqrt(dg_dr[0][3]*dg_dr[0][3] - dg_dr[0][0]*dg_dr[3][3])) / dg_dr[3][3];
  double angle = acos(sqrt(gxx(chi, r, theta, 1.0)) * sqrt(-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega) / (1.0 + b * Omega) * g[2][2] * vars[3]);


  // E = -(g[0][0] + Omega * g[0][3]) / sqrt(-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega);
  // Lz = (g[0][3] + Omega * g[3][3]) / sqrt(-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega);
  // kt = (E*g[3][3] + Lz*g[0][3]) / (g[0][3]*g[0][3] - g[0][0]*g[3][3]);

  // printf(" %f %f %f\n", (E*g[3][3] + Lz*g[0][3]) / (g[0][3]*g[0][3] - g[0][0]*g[3][3]), 1 / sqrt(-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega), kt1*E1);

  // printf(" r=%f, Thomas's=%f, Cosimo's=%f %e %e\n", r, acos(q / (r * kt)), angle, th1, th2);

  // printf(" %f %f, \n", q, g[2][2]*kth1 / kt1);

  // printf(" %f %f %f\n", r, acos(q / (r * kt)),  acos(q / (r * cal_kt(r, theta))));


  // printf(" q=%f kt=%f acos=%f", q, kt, q / (r * kt));
  // printf(" r = %f, q = %f, kt = %f, angle = %f\n", r, q, kt, acos(q / (r * kt)));

  // return acos(q / (r * kt));
  return angle;
}

// double cal_gi(double delta, double theta, double h, double r) {
//   double a=chi;
//   return cal_kt(r, theta)*sqrt(h*h-2*h+a*a)/sqrt(h*h + a*a);
// }

double cal_kt(double r, double theta) {
  double g_tt, g_pp, g_tp, denom, kt, a = chi, M=1.0;
  double r2, a2, s2, c2;
  r2 = r*r;
  a2 = a*a;
  s2 = sin(theta)*sin(theta);
  c2 = cos(theta)*cos(theta);
  g_tt = -1.0 + (2.0*M*r)/(r2 + a2*c2);
  // g_tt = -1. + (2.*M*r)/(pow(r,2) + pow(a,2)*pow(cos(theta),2));
  // g_pp = pow(sin(theta),2)*(pow(a,2) + pow(r,2) + (2*pow(a,2)*M*r*pow(sin(theta),2))/(pow(r,2) + pow(a,2)*pow(cos(vars[1]),2)));
  g_pp = s2*(a2 + r2 + (2.0*a2*M*r*s2)/(r2 + a2*c2));
  g_tp = (-2.0*a*M*r*s2)/(r2 + a2*c2);
  // g_tp = (-2.0*a*M*r*pow(sin(theta),2))/(pow(r,2) + pow(a,2)*pow(cos(theta),2));  

  
  denom = (g_tt*g_pp-g_tp*g_tp);
  kt = -(g_pp+b*g_tp)/denom;
  // printf(" gtt=%f, gpp=%f, gtp=%f ", g_tt, g_pp, g_tp);
  return kt;
}


void cache() {
  int i = 0;
  for (;i<N;i++)
    oldvars[i] = vars[i];
  oldktphi[0]=ktphi[0];
  oldktphi[1]=ktphi[1];
}

void rcache() {
  int i = 0;
  for (;i<N;i++)
    vars[i] = oldvars[i];
  ktphi[0] = oldktphi[0];
  ktphi[1] = oldktphi[1];
}

double R(double vars[5]) {
  double x, y, z;
  
  x = vars[0]*sin(vars[1])*cos(vars[4]);
  y = vars[0]*sin(vars[1])*sin(vars[4]);
  z = vars[0]*cos(vars[1]);

  return sqrt(x*x + y*y + z*z);
}

double interp_lin_1d(double ifac_r, double rlo, double rhi){
  return ifac_r*rhi + (1.0-ifac_r)*rlo;
}


double find_radii(double nvars[5], double ovars[5]) {
  // double x1, x2, y1, y2, z1, z2, x12, y12, z12, x, y;
  // x1 = ovars[0]*sin(ovars[1])*cos(ovars[4]);
  // y1 = ovars[0]*sin(ovars[1])*sin(ovars[4]);
  // z1 = ovars[0]*cos(ovars[1]);
  
  // x2 = nvars[0]*sin(nvars[1])*cos(nvars[4]);
  // y2 = nvars[0]*sin(nvars[1])*sin(nvars[4]);
  // z2 = nvars[0]*cos(nvars[1]);

  // // 0.5*sqrt(x1*x1 + y1*y1 + z1*z1) + 0.5*sqrt(x2*x2 + y2*y2 + z2*z2)

  // x12 = x1 - x2;
  // z12 = z1 - z2;
  // y12 = y1 - y2;

  // x = (x12 * z2 - z12 * x2) / z12;
  // y = (y12 * z2 - z12 * y2) / z12;
  // if (nvars[0] >= d)
  //   return d;

  return (ovars[0]+nvars[0])/2.0;//sqrt(x * x + y * y);

}
