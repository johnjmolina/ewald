#include "test_aux.h"
/*
  Hundred random point charges
 */
int main(int argc, char *argv[]){
  const char* name="q_100_rand";

  num = 100;
  init(num);
  set_cubic_box(10.0);

  double totalq = 0.0;
  double dmyq = -1.0;
  for(int i = 0; i < num; i++){
    q[i] = dmyq;
    dmyq*= -1.0;
    totalq += q[i];
  }
  assert(zero_mp(totalq));

  for(int i = 0; i < num; i++){
    r[i][0] = RAx(boxlen);
    r[i][1] = RAx(boxlen);
    r[i][2] = RAx(boxlen);
    dval[i] = 0.0;
    for(int d = 0; d < DIM; d++){
      mu[i][d] = 0.0;
    }
  }
  ndirect = 16;
  compute_all(true, true, name);
  free();
  delete ewald_sum;
  return 0;
}
