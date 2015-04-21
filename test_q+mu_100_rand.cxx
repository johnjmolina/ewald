#include "test_aux.h"
/*
  Hundred randomly oriented dipoles
 */
int main(int argc, char *argv[]){
  const char* name="q+mu_100_rand";

  num = 100;
  init(num);
  set_cubic_box(10.0);

  double qval = -1.0;
  double totalq = 0.0;
  for(int i = 0; i < num; i++){
    r[i][0] = RAx(boxlen);
    r[i][1] = RAx(boxlen);
    r[i][2] = RAx(boxlen);
    q[i]    = qval;
    dval[i] = 1.0;
    qval *= -1.0;
    totalq += q[i];
  }
  assert(zero_mp(totalq));
  for(int i = 0; i < num; i++){
    random_dipole(dval[i], mu[i]);
  }
  ndirect = 16;
  compute_all(true, true, false, name);
  free();
  delete ewald_sum;
  return 0;
}
