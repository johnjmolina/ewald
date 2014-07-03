#ifndef TEST_AUX_H
#define TEST_AUX_H
#include <time.h>
#include "ewald_gold.h"
#include "ewald.h"
#include "rigid_body.h"

int num;
double *dval;

double **r;
double *q;
double **mu;
double ***theta;
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

void init(const int &num){
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
  
  dval = (double*) malloc(num * sizeof(double));
  r = (double**) alloc_2d_double(num, DIM);
  q = (double*) alloc_1d_double(num);
  mu= (double**) alloc_2d_double(num, DIM);
  theta = (double***) alloc_3d_double(num, DIM, DIM);
  {
    for(int i = 0; i < num; i++){
      q[i] = dval[i] = 0.0;
    }
    double* rr = r[0];
    double* pp = mu[0];
    for(int i = 0; i < num*DIM; i++){
      rr[i] = pp[i] = 0.0;
    }
    double* tt = theta[0][0];
    for(int i = 0; i < num*DIM*DIM; i++){
      tt[i] = 0.0;
    }
  }

  force  = (double**) alloc_2d_double(num, DIM);
  torque = (double**) alloc_2d_double(num, DIM);
  efield = (double**) alloc_2d_double(num, DIM);
  force_gold  = (double**) alloc_2d_double(num, DIM);
  torque_gold = (double**) alloc_2d_double(num, DIM);
  efield_gold = (double**) alloc_2d_double(num, DIM);
}
void free(){
  free(dval);
  free_2d_double(r);
  free_1d_double(q);
  free_2d_double(mu);
  free_3d_double(theta);
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

void show_results(const int &n, const double &energy, 
                  double const* const* force, double const* const* torque, double const* const* efield,
                  FILE* fout){
  fprintf(fout, "%14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
          energy, 
	  force[0][0], force[0][1], force[0][2],
	  torque[0][0], torque[0][1], torque[0][2],
	  efield[0][0], efield[0][1], efield[0][2]
          );
  if(n > 1){
  fprintf(fout, "%14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
          energy, 
	  force[1][0], force[1][1], force[1][2],
	  torque[1][0], torque[1][1], torque[1][2],
	  efield[1][0], efield[1][1], efield[1][2]
          );
  }
}
bool load_gold(const char* save_buffer){
  char buffer[256];
  bool loaded=false;
  sprintf(buffer, "%s_gold.dat", save_buffer);
  FILE* fsave = fopen(buffer, "r");
  if(fsave != NULL){
    loaded = true;
    int dmy_n;
    fscanf(fsave, "%d", &dmy_n);
    assert(dmy_n == num);
    fscanf(fsave, "%lf", &energy_gold);
    for(int i = 0; i < num; i++){
      fscanf(fsave, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", 
             &force_gold[i][0], &force_gold[i][1], &force_gold[i][2],
             &torque_gold[i][0], &torque_gold[i][1], &torque_gold[i][2],
             &efield_gold[i][0], &efield_gold[i][1], &efield_gold[i][2]
             );
    }
  }
  return loaded;
}
void compute_gold(const char* save_buffer){
  if(!load_gold(save_buffer)){
    energy_gold = 0.0;
    ewald_direct_sum(energy_gold, force_gold, torque_gold, efield_gold, ndirect,
                     num, *cell, r, q, mu, stderr, save_buffer);
  }else{
    fprintf(stderr, "******* Reference results loaded from %s_gold.dat\n", save_buffer);
  }
  show_results(num, energy_gold, force_gold, torque_gold, efield_gold, stderr);
}
void compute_all(const bool& charge, const bool& dipole, const bool& quadrupole, const char* save_buffer){

  double* dmy_q     = (charge ? q : NULL);
  double* dmy_mu    = (dipole ? mu[0]: NULL);
  double* dmy_theta = (quadrupole ? theta[0][0]: NULL);

  fprintf(stderr, "******* Direct Sum Calculation\n");  
  compute_gold(save_buffer);

  fprintf(stderr, "\n****** Ewald Calculation\n");
  epsilon = 1.0;
  ewald_sum = new ewald(cell, alpha, epsilon, delta, conv, num,
                        charge, dipole, quadrupole);
  
  fprintf(stderr, "\t epsilon = 1 (vacuum)\n");
  epsilon = 1.0;
  ewald_sum -> reset_boundary(epsilon);
  ewald_sum -> compute(Ewald_energy, force[0], torque[0], efield[0], r[0], dmy_q, dmy_mu, dmy_theta);
  check_convergence(energy_gold, Ewald_energy[0], rmstol);
  show_results(num, Ewald_energy[0], force, torque, efield, stderr);
  
  fprintf(stderr, "\t epsilon = inf (tinfoil)\n");
  epsilon = -1.0;
  ewald_sum -> reset_boundary(epsilon);
  ewald_sum -> compute(Ewald_energy, force[0], torque[0], efield[0], r[0], dmy_q, dmy_mu, dmy_theta);
  check_convergence(energy_gold, Ewald_energy[0], rmstol2);
  show_results(num, Ewald_energy[0], force, torque, efield, stderr);
  
  print_convergence(rmstol, rmstol2, stderr);
}
#endif
