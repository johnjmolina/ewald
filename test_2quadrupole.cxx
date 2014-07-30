#include "test_aux.h"
/*
  Test on quadrupoles
 */
int main(int argc, char *argv[]){
  const char* name="2quadrupole";

  num = 2;
  init(num);
  set_cubic_box(10.0);
  
  for(int i = 0; i < num; i++){
    q[i]    = 0.0;
    dval[i] = 0.0;
    mu[i][0] = mu[i][1] = mu[i][2] = 0.0;
    theta[i][0][0] = theta[i][0][1] = theta[i][0][2] = 0.0;
    theta[i][1][0] = theta[i][1][1] = theta[i][1][2] = 0.0;
    theta[i][2][0] = theta[i][2][1] = theta[i][2][2] = 0.0;
  } 	
  r[0][0] = 4.5;
  r[0][1] = 5.0;
  r[0][2] = 5.0;

  r[1][0] = 5.5;
  r[1][1] = 5.0;
  r[1][2] = 5.0;

  double dmy_q = 1.0;
  double tot_q = 0.0;
  double tr_theta = 0.0;
  for(int i = 0; i < num; i++){
    tr_theta = random_quadrupole(1.0, theta[i][0]);
    for(int d = 0; d < DIM; d++){
      theta[i][d][d] -= (tr_theta / 3.0);
    }
    assert(zero_mp(theta[i][0][0] + theta[i][1][1] + theta[i][2][2]));
  }
  assert(tot_q == 0.0);

  compute_all(true, true, true, name);
  free();
  delete ewald_sum;
  return 0;
}
