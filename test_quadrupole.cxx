#include "test_aux.h"
/*
  Test on quadrupoles
 */
int main(int argc, char *argv[]){
  const char* name="24quadrupole";

  num = 24;
  init(num);
  set_cubic_box(10.0);
  
  for(int i = 0; i < num; i++){
    r[i][0] = RAx(10.0);
    r[i][1] = RAx(10.0);
    r[i][2] = RAx(10.0);

    q[i]    = 0.0;
    dval[i] = 0.0;
    mu[i][0] = mu[i][1] = mu[i][2] = 0.0;
    theta[i][0][0] = theta[i][0][1] = theta[i][0][2] = 0.0;
    theta[i][1][0] = theta[i][1][1] = theta[i][1][2] = 0.0;
    theta[i][2][0] = theta[i][2][1] = theta[i][2][2] = 0.0;
  } 	

  double dmy_q = 1.0;
  double tot_q = 0.0;
  double tr_theta = 0.0;
  for(int i = 0; i < num; i++){

    //charges
    q[i] = dmy_q;
    tot_q += q[i];
    dmy_q *= (-1.0);

    //dipoles
    dval[i] = 1.0;
    random_dipole(dval[i], mu[i]);

    //quadrupoles
    dval[i] = 1.0;
    tr_theta = random_quadrupole(dval[i], theta[i][0]);
  }
  assert(tot_q == 0.0);

  compute_all(true, true, true, name);
  free();
  delete ewald_sum;
  return 0;
}
