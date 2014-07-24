#include "test_aux.h"
/*
  Test on quadrupoles
 */
int main(int argc, char *argv[]){
  const char* name="10quadrupole";

  num = 2;
  init(num);
  set_cubic_box(10.0);
  
  r[0][0] = 4.5;
  r[1][0] = 5.5;
  r[0][1] = r[1][1] = 5.0;
  r[0][2] = r[1][2] = 5.0;
  for(int i = 0; i < num; i++){
    q[i]    = 0.0;
    mu[i][0] = mu[i][1] = mu[i][2] = 0.0;
  } 		
  for(int i = 0; i < num; i++){
    for(int d = 0; d < DIM; d++){
      for(int e = d; e < DIM; e++){
        theta[i][d][e] = theta[i][e][d] = RAx(1.0);
      }	
    }
    double tti = (theta[i][0][0] + theta[i][1][1] + theta[i][2][2]);
    fprintf(stderr, "tr %d : %.8E\n", i, tti);
    /*for(int d = 0; d < DIM; d++){
      theta[i][d][d] += (-tti/3.0);
      }*/
  }

  ndirect = 12;
  compute_all(false, false, true, name);
  free();
  delete ewald_sum;
  return 0;
}
