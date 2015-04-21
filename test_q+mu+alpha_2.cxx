#include "test_aux.h"

/*
  Two point charges + Two induced dipoles
 */
#include "test_aux.h"
void reset_dipole(double** mu_ind){
  for(int i = 0; i < num; i++){
    double* imu = mu_ind[i];
    imu[0] = imu[1] = imu[2] = 0.0;
  }
}
void induced_dipole(double** mu_ind, 
                    double const* const* field_ind){
  for(int i = 0; i < num; i++){
    double* imu = mu_ind[i];
    const double* iee = field_ind[i];
    const double* iaa = polar[i][0];
    imu[0] = iaa[0]*iee[0] + iaa[1]*iee[1] + iaa[2]*iee[2];
    imu[1] = iaa[3]*iee[0] + iaa[4]*iee[1] + iaa[5]*iee[2];
    imu[2] = iaa[6]*iee[0] + iaa[7]*iee[2] + iaa[8]*iee[2];
  }
}
void total_dipole(double** mu_ind,
                  double const* const* mu_perm){
  for(int i = 0; i < num; i++){
    double* imu = mu_ind[i];
    const double* imu0= mu_perm[i];

    imu[0] += imu0[0];
    imu[1] += imu0[1];
    imu[2] += imu0[2];
  }
}
void induction_energy(double &energy, double const* const* field_ind){
  double dmy = 0.0;
  for(int i = 0; i < num; i++){
    const double* iee = field_ind[i];
    const double* iaa = polar[i][0];
    dmy += iee[0]*(iaa[0]*iee[0] + iaa[1]*iee[1] + iaa[2]*iee[2])
      + iee[1]*(iaa[3]*iee[0] + iaa[4]*iee[1] + iaa[5]*iee[2])
      + iee[2]*(iaa[6]*iee[0] + iaa[7]*iee[1] + iaa[8]*iee[2]);
  }
  energy += dmy/2.0;
}
int main(int argc, char *argv[]){
  const char* name="q+mu+alpha_2";

  num=4;
  init(num);
  for(int i = 0; i < num; i++){
    q[i] = 0.0;
    mu0[i][0] = mu0[i][1] = mu0[i][2] = 0.0;
    mu[i][0] = mu[i][1] = mu[i][2] = 0.0;
    for(int d = 0; d < DIM; d++){
      theta[i][d][0] = theta[i][d][1] = theta[i][d][2] = 0.0;
      polar[i][d][0] = polar[i][d][1] = polar[i][d][2] = 0.0;
    }
  }

  {
    int id;
    id = 0;
    r[id][0] = 4.5;
    r[id][1] = 5.0;
    r[id][2] = 5.0;
    q[id]    = 1.0;
    
    id = 1;
    r[id][0] = 5.5;
    r[id][1] = 5.0;
    r[id][2] = 5.0;
    q[id]    = -1.0;
    
    id = 2;
    r[id][0] = 4.5;
    r[id][1] = 6.0;
    r[id][2] = 5.0;
    polar[id][0][0] = 0.1;
    
    id = 3;
    r[id][0] = 5.5;
    r[id][1] = 6.0;
    r[id][2] = 6.0;
    polar[id][0][0] = polar[id][1][1] = polar[id][2][2] = 0.1;
  }

  set_cubic_box(10.0);

  //gold calculation
  const int niter = 3;
  reset_dipole(mu);
  print_gold = false;
  for(int i = 0; i < niter; i++){
    energy_gold = 0.0;

    ewald_direct_sum(energy_gold, force_gold, torque_gold,
                     efield_gold, efield_grad_gold,
                     ndirect, num, *cell, r, q, mu, theta, stderr, name);
    induced_dipole(mu, efield_gold);
    total_dipole(mu, mu0);
    induction_energy(energy_gold, efield_gold);
    fprintf(stderr, "#ALPHA iter %d\n", i);
    show_results(num, energy_gold, force_gold, torque_gold, efield_gold,
                 efield_grad_gold, stderr);
    fprintf(stderr, "\n");
  }
  print_gold = true;
  energy_gold = 0.0;
  ewald_direct_sum(energy_gold, force_gold, torque_gold,
                   efield_gold, efield_grad_gold,
                   ndirect, num, *cell, r, q, mu, theta, stderr, name);
  induction_energy(energy_gold, efield_gold);
  show_results(num, energy_gold, force_gold, torque_gold, efield_gold,
               efield_grad_gold, stderr);
  fprintf(stderr, "induction =  %.8f %.8f %.8f %.8f %.8f %.8f\n", 
          mu[2][0], mu[2][1], mu[2][2],
          mu[3][0], mu[3][1], mu[3][2]);


  //ewald calculation
  reset_dipole(mu);
  epsilon = 1.0;
  ewald_sum = new ewald(cell, alpha, epsilon, delta, conv, num, true, true, false);
  Ewald_energy[0] = Ewald_energy[1] = Ewald_energy[2] = Ewald_energy[3] = Ewald_energy[4] = 0.0;
  for(int i = 0; i <= niter; i++){
    Ewald_energy[0] = Ewald_energy[1] = Ewald_energy[2] = Ewald_energy[3] = Ewald_energy[4] = 0.0;
    ewald_sum -> compute(Ewald_energy, force[0], torque[0], efield[0], efield_grad[0][0],
                         r[0], q, mu[0], theta[0][0], "alpha_ewald.dat");
    ewald_sum -> compute_mu_induced(mu[0], polar[0][0], efield[0]);
    ewald_sum -> compute_upol(Ewald_energy[0], polar[0][0], efield[0]);

    total_dipole(mu, mu0);
    fprintf(stderr, "#ALPHA iter %d\n", i);
    show_results(num, Ewald_energy[0], force, torque, efield, efield_grad, stderr);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "induction =  %.8f %.8f %.8f %.8f %.8f %.8f\n", 
          mu[2][0], mu[2][1], mu[2][2],
          mu[3][0], mu[3][1], mu[3][2]);
  check_convergence(energy_gold, Ewald_energy[0], rmstol);
  print_convergence(rmstol, rmstol2, stderr);

  free();
  delete ewald_sum;
  return 0;
}
