#include "test_aux.h"
int main(int argc, char* argv[]){
  num = 1;
  init_dipole(num);
  
  p[0].x[0] = 5.0;
  p[0].x[1] = 5.0;
  p[0].x[2] = 5.0;
  dval[0] = 1.0;
  qtn_init(p[0].q, 1.0, 0.0, 0.0, 0.0);
  double dmy_mu[DIM];
  rigid_body_rotation(dmy_mu, ex, p[0].q, BODY2SPACE);
  for(int d = 0; d < DIM; d++){
    r[0][d] = p[0].x[d];
    mu[0][d] = dval[0] * dmy_mu[d];
  }
  fprintf(stderr, "******* Direct Sum Calculation\n");
  energy_gold = 0.0;
  dipole_direct_sum(energy_gold, force_gold, torque_gold, efield_gold, ndirect,
                    num, *cell, r, mu, stderr);
  fprintf(stderr, "%14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
          energy_gold, 
	  force_gold[0][0], force_gold[0][1], force_gold[0][2],
	  torque_gold[0][0], torque_gold[0][1], torque_gold[0][2],
	  efield_gold[0][0], efield_gold[0][1], efield_gold[0][2]
          );


  fprintf(stderr, "\n******* Ewald Calculation\n");
  epsilon = 1.0;
  ewald_sum = new ewald(cell, alpha, epsilon, delta, conv, num,
                               false, true, false);

  fprintf(stderr, "\t epsilon = 1 (vacuum)\n");
  ewald_sum -> reset_boundary(epsilon);
  ewald_sum -> compute(Ewald_energy, force[0], torque[0], efield[0], r[0], NULL, mu[0], NULL);
  check_convergence(energy_gold, Ewald_energy[0], rmstol);
  fprintf(stderr, "%14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
          Ewald_energy[0], force[0][0], force[0][1], force[0][2],
	  torque[0][0], torque[0][1], torque[0][2],
	  efield[0][0], efield[0][1], efield[0][2]
          );

  fprintf(stderr, "\t epsilon = inf (tinfoil)\n");
  epsilon = -1.0;
  ewald_sum -> reset_boundary(epsilon);
  ewald_sum -> compute(Ewald_energy, force[0], torque[0], efield[0], r[0], NULL, mu[0], NULL);
  check_convergence(energy_gold, Ewald_energy[0], rmstol2);
  fprintf(stderr, "%14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
          Ewald_energy[0], force[0][0], force[0][1], force[0][2],
	  torque[0][0], torque[0][1], torque[0][2],
	  efield[0][0], efield[0][1], efield[0][2]
          );

  print_convergence(rmstol, rmstol2, stderr);
  delete ewald_sum;
  free_dipole();

  return 0;
}
