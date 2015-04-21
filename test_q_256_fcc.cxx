#include "test_aux.h"
void fccLattice(int ncells)
{
  int Ix,Iy,Iz,i,Ip,In;

  i=0;
  Ip=0;
  In=0;

  for(Ix=0;Ix<ncells;Ix++){
    for(Iy=0;Iy<ncells;Iy++){
      for(Iz=0;Iz<ncells;Iz++){
        i=Iz+Iy*ncells+Ix*ncells*ncells;

        r[i][0]=Ix+0.5;
        r[i][1]=Iy+0.5;
        r[i][2]=Iz+0.5;

        if(((Ix+Iy+Iz)%2)==0){
          q[i]=1.0;
          Ip++;
        }else{
          q[i]=-1.0;
          In++;
        }
	mu[i][0] = mu[i][1] = mu[i][2] = 0.0;
      }
    }
  }
}
/*
  Hundred random point charges
 */
int main(int argc, char *argv[]){
  const char* name="256fcc";
  const int ncells = 8;
  const int num_fcc = ncells*ncells*ncells;
  num = num_fcc;
  init(num);
  set_cubic_box(float(ncells));
  fccLattice(ncells);
  ndirect = 12;
  compute_all(true, true, true, name);
  free();
  delete ewald_sum;
  return 0;
}
