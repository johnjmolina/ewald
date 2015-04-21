#include "test_aux.h"
/*
  Two point charges
 */
int main(int argc, char *argv[]){
  const char* name="q_2";

  num = 2;
  init(num);

  r[0][0] = 4.5;
  r[1][0] = 5.5;
  r[0][1] = r[1][1] = 5.0;
  r[0][2] = r[1][2] = 5.0;
  q[0] = 1.0;
  q[1] = -1.0;
  dval[0] = dval[1] = 0.0;
  mu[0][0] = mu[0][1] = mu[0][2] = 0.0;
  mu[1][0] = mu[1][1] = mu[1][2] = 0.0;
  set_cubic_box(10.0);
  compute_all(true, true, false, name);
  free();
  delete ewald_sum;
  return 0;
}
