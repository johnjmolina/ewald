#include "test_aux.h"
/*
  Two randomly oriented dipoles
 */
int main(int argc, char *argv[]){
  const char* name="2dipole_random";

  num = 2;
  init(num);

  r[0][0] = r[0][1] = r[0][2] = 4.5;
  r[1][0] = r[1][1] = r[1][2] = 5.5;
  q[0] = q[1] = 0.0;
  dval[0] = dval[1] = 1.0;

  random_dipole(dval[0], mu[0]);
  random_dipole(dval[1], mu[1]);

  set_cubic_box(10.0);
  compute_all(true, true, false, name);
  free();
  delete ewald_sum;
  return 0;
}
