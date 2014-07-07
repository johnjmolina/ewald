#include "ewald_gold.h"

inline void pair_interaction(const double rij[DIM],
                             const double &qi, const double &qj,
                             const double mui[DIM], const double muj[DIM],
                             double &energy, 
                             double force[DIM], double torque[DIM], double field[DIM]
                             ){
  energy = 0.0;
  force[0] = force[1] = force[2] = 0.0;
  torque[0] = torque[1] = torque[2] = 0.0;
  field[0] = field[1] = field[2] = 0.0;

  double rr = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
  double dri = 1.0/rr;
  double dr2i = dri*dri;
  
  double Br = dr2i*dri;
  double Cr = 3.0*Br*dr2i;
  double Dr = 5.0*Cr*dr2i;
  
  { //charge
    energy += qi*qj*dri;
    for(int d = 0; d < DIM; d++){
      force[d] += qi*qj*Br*rij[d];
      field[d] += qj*Br*rij[d];
    }
  }
  {  //dipole
    double dot_mui_r = v_inner_prod(mui, rij);
    double dot_muj_r = v_inner_prod(muj, rij);
    double dot_mui_muj = v_inner_prod(mui, muj);
    double cross_mui_r[DIM], cross_muj_r[DIM], cross_mui_muj[DIM];
    v_cross(cross_mui_r, mui, rij);
    v_cross(cross_muj_r, muj, rij);
    v_cross(cross_mui_muj, mui, muj);
    
    //charge-dipole
    //    energy += (qi*dot_muj_r - qj*dot_mui_r)*Br;

    //dipole-dipole
    energy += (Br*dot_mui_muj - Cr*dot_mui_r*dot_muj_r);
    for(int d = 0; d < DIM; d++){
      force[d]  += (Cr*(dot_mui_muj*rij[d] + dot_muj_r*mui[d] + dot_mui_r*muj[d])
                    - Dr*dot_mui_r*dot_muj_r*rij[d]);
      torque[d] += (-Br*cross_mui_muj[d] + Cr*cross_mui_r[d]*dot_muj_r);
      field[d]  += (-Br*muj[d] + Cr*rij[d]*dot_muj_r);
    }
  }
}
void ewald_minimum_image(double & energy,
                         double **force,
                         double **torque,
                         double **field,
                         const int &Nparticles,
                         const parallelepiped &rcell,
                         double const* const* r,
                         double const* q,
                         double const* const* mu,
                         FILE* fout,
                         const char* save_buffer
                         ){
  double dmy_energy;
  double dmy_force[DIM], dmy_torque[DIM], dmy_field[DIM];
  double rij_mi[DIM];
  for(int i = 0; i < Nparticles; i++){
    const double qi   = q[i];
    const double* mui = mu[i];

    for(int j = 0; j < Nparticles; j++){
      const double qj   = q[j];
      const double* muj = mu[j];

      if(i != j){
        rcell.distance_MI(r[i], r[j], rij_mi);
        pair_interaction(rij_mi, qi, qj, mui, muj, 
                         dmy_energy, dmy_force, dmy_torque, dmy_field);
        
        energy += (dmy_energy/2.0);
        for(int d = 0; d < DIM; d++){
          force[i][d]  += dmy_force[d];
          torque[i][d] += dmy_torque[d];
          field[i][d]  += dmy_field[d];
        }
      }
    }
  }

  {
    char buffer[256];
    sprintf(buffer, "%s_mi.dat", save_buffer);
    FILE* fsave = filecheckopen(buffer, "w");
    fprintf(fsave, "%d\n", Nparticles);
    fprintf(fsave, "%20.12E\n", energy);
    for(int i = 0; i < Nparticles; i++){
      fprintf(fsave, "%20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E\n",
              force[i][0], force[i][1], force[i][2],
              torque[i][0], torque[i][1], torque[i][2],
              field[i][0], field[i][1], field[i][2]
              );
    }
    fclose(fsave);
  }
}
void ewald_direct_sum(double &energy, 
                      double **force, 
                      double **torque, 
                      double **field, 
                      const int &ncut,
                      const int &Nparticles,
                      const parallelepiped &rcell,
                      double const* const* r,
                      double const* q,
                      double const* const* mu,
                      FILE* fout,
                      const char* save_buffer
                      ){
  vector<double> nmag;
  vector< vector< vector<int> > > nlist;
  long int ntot;
  int lmax, mmax, nmax;
  int icell[DIM];
  int &ll = icell[0];
  int &mm = icell[1];
  int &nn = icell[2];

  double** shell_force = (double**)alloc_2d_double(Nparticles, DIM);
  double** shell_torque = (double**)alloc_2d_double(Nparticles, DIM);
  double** shell_field = (double**)alloc_2d_double(Nparticles, DIM);
  double shell_energy;
  double dmy_force[DIM];
  double dmy_torque[DIM];
  double dmy_field[DIM];
  double lbox[DIM];
  double dmy_energy;
  double rij_mi[DIM], rij[DIM];
  clock_t start_t, end_t;
  double cpu_t, cpu_t2;


  // generate spherical shells of periodic cells 
  rcell.get_lengths(lbox[0], lbox[1], lbox[2]);
  nshell_list(ncut, lbox[0], lbox[1], lbox[2], 
	      lmax, mmax, nmax, ntot, nmag, nlist, fout);
  start_t = clock();
  ewald_direct_sum_naive(energy, force, torque, field, lmax, mmax, nmax, Nparticles, rcell, 
                         r, q, mu, fout, save_buffer);
  end_t = clock();
  cpu_t2 = ((double)end_t - start_t)/CLOCKS_PER_SEC;

  start_t = clock();
  energy = shell_energy = 0.0;
  for(int i = 0; i < Nparticles; i++){
    for(int d = 0; d < DIM; d++){
      force[i][d] = shell_force[i][d] = 0.0;
      torque[i][d] = shell_torque[i][d] = 0.0;
      field[i][d] = shell_field[i][d] = 0.0;
    }
  }

  for(int nshell = 0; nshell < (int)nlist.size(); nshell++){
    vector< vector<int> > &shell = nlist.at(nshell);

    for(int ncell = 0; ncell < shell.size(); ncell++){
      vector<int> &cell = shell.at(ncell);
      icell[0] = cell.at(0);
      icell[1] = cell.at(1);
      icell[2] = cell.at(2);
      
      for(int i = 0; i < Nparticles; i++){
        const double qi = q[i];
	const double* mui = mu[i];
	
	for(int j = 0; j < Nparticles; j++){
          const double qj = q[j];
	  const double* muj = mu[j];
	  if(!(i == j && ll == 0 && mm == 0 && nn == 0)){
            
	    rcell.distance_MI(r[i], r[j], rij_mi);
	    for(int d = 0; d < DIM; d++){
	      rij[d] = rij_mi[d] + (double)icell[d]*lbox[d];
	    }
            pair_interaction(rij, qi, qj, mui, muj, dmy_energy, dmy_force, dmy_torque, dmy_field);
	    
	    shell_energy += (dmy_energy/2.0);
	    for(int d = 0; d < DIM; d++){
	      shell_force[i][d] += dmy_force[d];
	      shell_torque[i][d] += dmy_torque[d];
	      shell_field[i][d] += dmy_field[d];
	    }
	    
	  }//valid pairs
	  
	}//j
      }//i
      
    }//cells in nshell

    //update cumulative values
    energy += shell_energy;
    for(int i = 0; i < Nparticles; i++){
      for(int d = 0; d < DIM; d++){
	force[i][d] += shell_force[i][d];
	torque[i][d] += shell_torque[i][d];
	field[i][d] += shell_field[i][d];
      }
      if(i == 0){
	fprintf(fout, 
		"%3d %10d %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
		nshell, (int) shell.size(),
		shell_energy, energy,
		force[0][0], force[0][1], force[0][2],
		torque[0][0], torque[0][1], torque[0][2],
		field[0][0], field[0][1], field[0][2]);
      }
      shell_energy = 0.0;
      for(int d = 0; d < DIM; d++){
	shell_force[i][d] = 0.0;
	shell_torque[i][d] = 0.0;
	shell_field[i][d] = 0.0;
      }
    }//i

  }//shells

  {
    char buffer[256];
    sprintf(buffer, "%s_gold.dat", save_buffer);
    FILE* fsave = filecheckopen(buffer, "w");
    fprintf(fsave, "%d\n", Nparticles);
    fprintf(fsave, "%20.12E\n", energy);
    for(int i = 0; i < Nparticles; i++){
      fprintf(fsave, "%20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E\n", 
              force[i][0], force[i][1], force[i][2],
              torque[i][0], torque[i][1], torque[i][2],
              field[i][0], field[i][1], field[i][2]
              );
    }
    fclose(fsave);
  }
  end_t = clock();
  cpu_t = ((double)end_t - start_t)/CLOCKS_PER_SEC;
  fprintf(fout, "\tExecution Time: spherical (naive)= %12.5f (%12.5f)\n", cpu_t, cpu_t2);
}

void ewald_direct_sum_naive(double &energy, 
                            double **force, 
                            double **torque, 
                            double **field,
                            const int &lmax, 
                            const int &mmax, 
                            const int &nmax,
                            const int &Nparticles,
                            const parallelepiped &rcell,
                            double const* const* r,
                            double const* q,
                            double const* const* mu,
                            FILE* fout,
                            const char* save_buffer
                            ){
  int icell[DIM];
  double lbox[DIM];
  double dmy_force[DIM];
  double dmy_torque[DIM];
  double dmy_field[DIM];
  double dmy_energy;
  double rij_mi[DIM], rij[DIM];

  energy = 0.0;
  for(int i = 0; i < Nparticles; i++){
    for(int d = 0; d < DIM; d++){
      force[i][d] = 0.0;
      torque[i][d] = 0.0;
      field[i][d] = 0.0;
    }
  }
  rcell.get_lengths(lbox[0], lbox[1], lbox[2]);

  for(int ll = -lmax; ll <= lmax; ll++){
    for(int mm = -mmax; mm <= mmax; mm++){
      for(int nn = -nmax; nn <= nmax; nn++){
	icell[0] = ll;
	icell[1] = mm;
	icell[2] = nn;

	for(int i = 0; i < Nparticles; i++){
          const double qi   = q[i];
	  const double* mui = mu[i];

	  for(int j = 0; j < Nparticles; j++){
            const double qj   = q[j];
	    const double* muj = mu[j];
            
	    if(!(i == j && ll == 0 && mm == 0 && nn == 0)){
	      rcell.distance_MI(r[i], r[j], rij_mi);
	      for(int d = 0; d < DIM; d++){
		rij[d] = rij_mi[d] + (double)icell[d]*lbox[d];
	      }
              pair_interaction(rij, qi, qj, mui, muj, dmy_energy, dmy_force, dmy_torque, dmy_field);

	      energy+= (dmy_energy/2.0);
	      for(int d = 0; d < DIM; d++){
		force[i][d] += dmy_force[d];
		torque[i][d] += dmy_torque[d];
		field[i][d] += dmy_field[d];
	      }

	    }//valid pairs
	  }//j
	}//i

      }//nn
    }//mm
  }//ll
  {
    char buffer[256];
    sprintf(buffer, "%s_gold0.dat", save_buffer);
    FILE* fsave = filecheckopen(buffer, "w");
    fprintf(fsave, "%d\n", Nparticles);
    fprintf(fsave, "%20.12E\n", energy);
    for(int i = 0; i < Nparticles; i++){
      fprintf(fsave, "%20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E\n", 
              force[i][0], force[i][1], force[i][2],
              torque[i][0], torque[i][1], torque[i][2],
              field[i][0], field[i][1], field[i][2]
              );
    }
    fclose(fsave);
  }
}
