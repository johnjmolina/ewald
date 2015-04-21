#include "test_aux.h"
/*
  One dipole
 */
int main(int argc, char* argv[]){
  const char* name="mu_1";

  num = 1;
  init(num);
  
  r[0][0] = 5.0;
  r[0][1] = 5.0;
  r[0][2] = 5.0;
  q[0]    = 0.0;
  dval[0] = 1.0;
  mu[0][0] = 1.0;
  mu[0][1] = mu[0][2] = 0.0;
  set_cubic_box(10.0);
  compute_all(true, true, false, name);
  free();
  delete ewald_sum;
  return 0;
}
