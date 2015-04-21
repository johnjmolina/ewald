#include "test_aux.h"


const double mass_o    = 15.9994;
const double mass_h    =  1.00794;
const double bond_h2o  =  0.9572;
const double charge_o  = -0.8010;
const double charge_h  =  0.4005;
const double angle_h2o =  104.52*M_PI/180.0;
const double mass_h2o  =  mass_o + 2.0*mass_h;

const double z_h = (mass_o/mass_h2o)*bond_h2o*cos(angle_h2o/2.0);
const double y_h = 0.0;
const double x_h = bond_h2o*sin(angle_h2o/2.0);

const double z_o = -(2.0*mass_h/mass_h2o)*bond_h2o*cos(angle_h2o/2.0);
const double y_o = 0.0;
const double x_o = 0.0;

/*
  Two randomly oriented dipoles
 */
int main(int argc, char *argv[]){
  const char* name="q+mu+alpha_2_water";

  const int nmol = 2;
  const int nat  = 3;

  num = nat*nmol;
  init(num);

  double **rcm = alloc_2d_double(nmol, DIM);

  //molecule 1
  rcm[0][0] = rcm[0][1] = rcm[0][2] = -1.85;
  //molecule 1
  rcm[1][0] = rcm[1][1] = rcm[1][2] = 0.0;

  for(int i = 0; i < nmol; i++){
    //o
    int io  = i*nat;
    int ih1 = io + 1;
    int ih2 = io + 2;
    group_id[io] = group_id[ih1] = group_id[ih2] = i;

    q[io] = charge_o;    
    r[io][0] = rcm[i][0] + x_o;
    r[io][1] = rcm[i][1] + y_o;
    r[io][2] = rcm[i][2] + z_o;

    //h1
    q[ih1]   = charge_h;
    r[ih1][0] = rcm[i][0] + x_h;
    r[ih1][1] = rcm[i][1] + y_h;
    r[ih1][2] = rcm[i][2] + z_h;

    //h1
    q[ih2]    = charge_h;    
    r[ih2][0] = rcm[i][0] - x_h;
    r[ih2][1] = rcm[i][1] + y_h;
    r[ih2][2] = rcm[i][2] + z_h;

    

    fprintf(stderr, "MOLECULE %d:\n", i);
    fprintf(stderr, "O  : id=%2d, r=(%.5f %.5f %.5f), q=%.3f \n", group_id[io],
	    r[i*nat][0], r[i*nat][1], r[i*nat][2], q[io]);
    fprintf(stderr, "H1 : id=%2d, r=(%.5f %.5f %.5f), q=%.3f \n", group_id[ih1],
	    r[i*nat+1][0], r[i*nat+1][1], r[i*nat+1][2], q[ih1]);
    fprintf(stderr, "H2 : id=%2d, r=(%.5f %.5f %.5f), q=%.3f \n", group_id[ih2],
	    r[i*nat+2][0], r[i*nat+2][1], r[i*nat+2][2], q[ih2]);

    dval[i] = 0.0;
    mu[i][0] = mu[i][1] = mu[i][2] = 0.0;
  }

  set_cubic_box(10.0);
  compute_all(true, false, name);
  free();
  delete ewald_sum;
  return 0;
}
