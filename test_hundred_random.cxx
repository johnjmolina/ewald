#include "test_aux.h"
/*
  Hundred randomly oriented dipoles
 */
int main(int argc, char *argv[]){
  const char* name="100dipole_random";

  num = 100;
  init(num);

  for(int i = 0; i < num; i++){
    r[i][0] = RAx(boxlen);
    r[i][1] = RAx(boxlen);
    r[i][2] = RAx(boxlen);
    q[i]    = 0.0;
    dval[i] = 1.0;

    quaternion rQ;
    double dmy_mu[DIM];
    random_rqtn(rQ);
    rigid_body_rotation(dmy_mu, ez, rQ, BODY2SPACE);
    for(int d = 0; d < DIM; d++){
      mu[i][d] = dval[i] * dmy_mu[d];
    }
  }
  compute_all(false, true, false, name);
  free();
  delete ewald_sum;
  return 0;
}
