#include "test_aux.h"
/*
  Two parallel dipoles
 */
int main(int argc, char *argv[]){
  const char* name="2dipole";

  num = 2;
  init(num);

  r[0][0] = 4.5;
  r[1][0] = 5.5;
  r[0][1] = r[1][1] = 5.0;
  r[0][2] = r[1][2] = 5.0;
  q[0] = q[1] = 0.0;
  dval[0] = dval[1] = 1.0;

  quaternion rQ;
  double dmy_mu[DIM];
  for(int i = 0; i < num; i++){
    qtn_init(rQ, 1.0, 0.0, 0.0, 0.0);
    rigid_body_rotation(dmy_mu, ex, rQ, BODY2SPACE);
    for(int d = 0; d < DIM; d++){
      mu[i][d] = dval[i] * dmy_mu[d];
    }
  }
  compute_all(false, true, false, name);
  free();
  delete ewald_sum;
  return 0;
}
