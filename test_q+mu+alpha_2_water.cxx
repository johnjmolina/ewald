#include "test_aux.h"


const double mass_o    = 15.9994;
const double mass_h    =  1.00794;
const double bond_h2o  =  0.9572;
const double angle_h2o =  104.52*M_PI/180.0;
const double mass_h2o  =  mass_o + 2.0*mass_h;

const double z_h = (mass_o/mass_h2o)*bond_h2o*cos(angle_h2o/2.0);
const double y_h = 0.0;
const double x_h = bond_h2o*sin(angle_h2o/2.0);

const double z_o = -(2.0*mass_h/mass_h2o)*bond_h2o*cos(angle_h2o/2.0);
const double y_o = 0.0;
const double x_o = 0.0;

const double charge_o  = -0.8010;
const double charge_h  =  0.4005;

const double mu_o[DIM]  = {0.0, 0.0, -0.2455};
const double mu_h1[DIM] = {0.1482, 0.0, 0.0797};
const double mu_h2[DIM] = {-0.1482, 0.0, 0.0797};

const double alpha_o[DIM][DIM] = {{4.5315, 0.0000, 0.0000},
				  {0.0000, 4.1388, 0.0000},
				  {0.0000, 0.0000, 4.7255}};

const double alpha_h1[DIM][DIM] = {{1.4837, 0.0000, 0.6782},
				   {0.0000, 0.5099, 0.0000},
				   {0.6782, 0.0000, 0.9721}};

const double alpha_h2[DIM][DIM] = {{1.4837, 0.0000, -0.6782},
				   {0.0000, 0.5099, 0.0000},
				   {-0.6782, 0.0000, 0.9721}};

/*
  Two randomly oriented dipoles
 */
int main(int argc, char *argv[]){
  const char* name="q+mu+alpha_2_water";

  const int nmol = 2;
  const int nat  = 3;

  num = nat*nmol;
  init(num);

  double **rcm = alloc_2d_double(nmol, DIM);

  //molecule 1
  rcm[0][0] = rcm[0][1] = rcm[0][2] = -1.85;
  //molecule 1
  rcm[1][0] = rcm[1][1] = rcm[1][2] = 0.0;

  for(int i = 0; i < nmol; i++){
    //o
    int io  = i*nat;
    int ih1 = io + 1;
    int ih2 = io + 2;
    group_id[io] = group_id[ih1] = group_id[ih2] = i;

    q[io] = charge_o;    
    r[io][0] = rcm[i][0] + x_o;
    r[io][1] = rcm[i][1] + y_o;
    r[io][2] = rcm[i][2] + z_o;
    for(int d = 0; d < DIM; d++) mu0[io][d] = mu_o[d];
    for(int l = 0; l < DIM; l++) for(int m = 0; m < DIM; m++) polar[io][l][m] = alpha_o[l][m];

    //h1
    q[ih1]   = charge_h;
    r[ih1][0] = rcm[i][0] + x_h;
    r[ih1][1] = rcm[i][1] + y_h;
    r[ih1][2] = rcm[i][2] + z_h;
    for(int d = 0; d < DIM; d++) mu0[ih1][d] = mu_h1[d];
    for(int l = 0; l < DIM; l++) for(int m = 0; m < DIM; m++) polar[ih1][l][m] = alpha_h1[l][m];

    //h2
    q[ih2]    = charge_h;    
    r[ih2][0] = rcm[i][0] - x_h;
    r[ih2][1] = rcm[i][1] + y_h;
    r[ih2][2] = rcm[i][2] + z_h;
    for(int d = 0; d < DIM; d++) mu0[ih2][d] = mu_h2[d];
    for(int l = 0; l < DIM; l++) for(int m = 0; m < DIM; m++) polar[ih2][l][m] = alpha_h2[l][m];

    fprintf(stderr, "MOLECULE %d:\n", i);
    fprintf(stderr, "O  : id=%2d, r=(%.5f %.5f %.5f), q=%.3f \n", group_id[io],
	    r[i*nat][0], r[i*nat][1], r[i*nat][2], q[io]);
    fprintf(stderr, "H1 : id=%2d, r=(%.5f %.5f %.5f), q=%.3f \n", group_id[ih1],
	    r[i*nat+1][0], r[i*nat+1][1], r[i*nat+1][2], q[ih1]);
    fprintf(stderr, "H2 : id=%2d, r=(%.5f %.5f %.5f), q=%.3f \n", group_id[ih2],
	    r[i*nat+2][0], r[i*nat+2][1], r[i*nat+2][2], q[ih2]);

  }

  set_cubic_box(10.0);

  // gold calculation
  const int niter = 3;
  reset_dipole(mu);
  total_dipole(mu, mu0);
  print_gold = false;
  for(int i = 0; i < niter; i++){
    energy_gold = 0.0;
    ewald_direct_sum(energy_gold, force_gold, torque_gold,
		     efield_gold, efield_grad_gold,
		     ndirect, num, *cell, group_id, r, q, mu, stderr, name);

    induced_dipole(mu, efield_gold, polar);
    total_dipole(mu, mu0);
    induction_energy(energy_gold, efield_gold, polar);
    fprintf(stderr, "#ALPHA iter %d\n", i);
    show_results(num, energy_gold, force_gold, torque_gold, efield_gold,
		 efield_grad_gold, stderr);
    fprintf(stderr, "\n");
  }

  print_gold = true;
  energy_gold= 0.0;
  ewald_direct_sum(energy_gold, force_gold, torque_gold,
		   efield_gold, efield_grad_gold,
		   ndirect, num, *cell, group_id, r, q, mu, stderr, name);
  induction_energy(energy_gold, efield_gold, polar);
  show_results(num, energy_gold, force_gold, torque_gold, efield_gold,
	       efield_grad_gold, stderr);
  fprintf(stderr, "induction = %.8f %.8f %.8f\n",
	  mu[0][0], mu[0][1], mu[0][2]
	  );


  //ewald calculation
  reset_dipole(mu);
  total_dipole(mu, mu0);
  epsilon = 1.0;

  ewald_sum = new ewald(cell, alpha, epsilon, delta, conv, num, true, true);
  ewald_sum -> define_groups(group_id);
  for(int i = 0; i <= niter; i++){
    Ewald_energy[0] = Ewald_energy[1] = Ewald_energy[2] = Ewald_energy[3] = Ewald_energy[4] = 0.0;
    ewald_sum -> compute(Ewald_energy, force[0], torque[0], efield[0], efield_grad[0][0],
			 r[0], q, mu[0], "alpha_ewald.dat");

    ewald_sum -> compute_mu_induced(mu[0], polar[0][0], efield[0]);
    ewald_sum -> compute_upol(Ewald_energy[0], polar[0][0], efield[0]);

    total_dipole(mu, mu0);
    fprintf(stderr, "#ALPHA iter %d\n", i);
    show_results(num, Ewald_energy[0], force, torque, efield, efield_grad, stderr);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "induction = %.8f %.8f %.8f\n",
	  mu[0][0], mu[0][1], mu[0][2]
	  );
  check_convergence(energy_gold, Ewald_energy[0], rmstol);
  print_convergence(rmstol, rmstol2, stderr);
  
  
  free();
  delete ewald_sum;
  return 0;
}
