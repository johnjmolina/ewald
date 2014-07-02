#ifndef TEST_AUX_H
#define TEST_AUX_H
#include <time.h>
#include "variable.h"
#include "dipole_gold.h"
#include "ewald.h"
#include "rigid_body.h"

int num;
Particle *p;
double *dval;

double **r;
double **mu;
double **force;
double **torque;
double **efield;
double Ewald_energy[5];

int ndirect;
double energy_gold;
double **force_gold;
double **torque_gold;
double **efield_gold;
const double ex[DIM] = {1.0, 0.0, 0.0};
const double ey[DIM] = {0.0, 1.0, 0.0};
const double ez[DIM] = {0.0, 0.0, 1.0};

double boxlen;
double alpha;
double delta;
double conv;
double epsilon;

double rmstol[4];
double rmstol2[4];
double a[DIM], b[DIM], c[DIM];
parallelepiped *cell;
ewald *ewald_sum;

void init_dipole(const int &num){
  boxlen = 10.0;
  alpha  = 8.0;
  delta  = 1.0e-16;
  conv   = 0.51;
  epsilon= 1.0;
  ndirect= 24;
  
  a[0] = b[1] = c[2] = boxlen;
  a[1] = b[2] = c[0] = 0.0;
  a[2] = b[0] = c[1] = 0.0;
  cell = new parallelepiped(a, b, c);
  
  p = (Particle*) malloc(num * sizeof(Particle));
  dval = (double*) malloc(num * sizeof(double));
  r = (double**) alloc_2d_double(num, DIM);
  mu= (double**) alloc_2d_double(num, DIM);
  force  = (double**) alloc_2d_double(num, DIM);
  torque = (double**) alloc_2d_double(num, DIM);
  efield = (double**) alloc_2d_double(num, DIM);
  force_gold  = (double**) alloc_2d_double(num, DIM);
  torque_gold = (double**) alloc_2d_double(num, DIM);
  efield_gold = (double**) alloc_2d_double(num, DIM);
}
void free_dipole(){
  free(p);
  free(dval);
  free_2d_double(r);
  free_2d_double(mu);
  free_2d_double(force);
  free_2d_double(torque);
  free_2d_double(efield);
  free_2d_double(force_gold);
  free_2d_double(torque_gold);
  free_2d_double(efield_gold);
}
void check_convergence(const double ener_gold, const double &ener, double rms[4]){
  double dmy_f[DIM];
  double dmy_t[DIM];
  double dmy_e[DIM];

  for(int i = 0; i < DIM; i++){
    dmy_f[i] = dmy_t[i] = dmy_e[i] = 0.0;
  }
  rms[0] = rms[1] = rms[2] = rms[3] = 0.0;
  for(int i = 0; i < num; i++){
    v_add(dmy_f, force_gold[i], force[i], -1.0);
    v_add(dmy_t, torque_gold[i], torque[i], -1.0);
    v_add(dmy_e, efield_gold[i], efield[i], -1.0);

    rms[1] += v_sqnorm(dmy_f);
    rms[2] += v_sqnorm(dmy_t);
    rms[3] += v_sqnorm(dmy_e);
  }

  rms[0] = sqrt((ener_gold - ener)*(ener_gold-ener));
  rms[1] = sqrt(rms[1] / (double)num);
  rms[2] = sqrt(rms[2] / (double)num);
  rms[3] = sqrt(rms[3] / (double)num);
}
void print_convergence(const double rms[4], const double rms2[4], FILE* fout){
  fprintf(fout, "\n RMS convergence parameters: vacuum (conducting)\n");
  fprintf(fout, "Energy = %1.0E (%1.0E)\n", rms[0], rms2[0]);
  fprintf(fout, "Force  = %1.0E (%1.0E)\n", rms[1], rms2[1]);
  fprintf(fout, "Torque = %1.0E (%1.0E)\n", rms[2], rms2[2]);
  fprintf(fout, "Field  = %1.0E (%1.0E)\n", rms[3], rms2[3]);
}
#endif
