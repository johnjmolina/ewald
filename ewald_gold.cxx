#include "ewald_gold.h"

bool print_gold;
// Compute pair interaction energy, as well as
// force, torque, and efield on particle i due to particle j
inline void pair_interaction(const double rij[DIM],
                             const double &qi, const double &qj,
                             const double mui[DIM], const double muj[DIM],
                             double &energy, 
                             double force[DIM], 
			     double torque[DIM], 
			     double field[DIM],
			     double field_grad[DIM][DIM]
                             ){
  energy = 0.0;
  force[0] = force[1] = force[2] = 0.0;
  torque[0] = torque[1] = torque[2] = 0.0;
  field[0] = field[1] = field[2] = 0.0;
  field_grad[0][0] = field_grad[0][1] = field_grad[0][2] = 0.0;
  field_grad[1][0] = field_grad[1][1] = field_grad[1][2] = 0.0;
  field_grad[2][0] = field_grad[2][1] = field_grad[2][2] = 0.0;

  double rr = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
  double dri = 1.0/rr;
  double dr2i = dri*dri;
  
  double Br = dr2i*dri;
  double Cr = 3.0*Br*dr2i;
  double Dr = 5.0*Cr*dr2i;
  double Er = 7.0*Dr*dr2i;
  double Fr = 9.0*Er*dr2i;

  double mui_r = v_inner_prod(mui, rij);
  double muj_r = v_inner_prod(muj, rij);
  double mui_muj = v_inner_prod(mui, muj);

  { //charge interactions
    energy += qi*qj*dri;
    for(int d = 0; d < DIM; d++){
      force[d]  += qi*qj*Br*rij[d];
      field[d]  += qj*Br*rij[d];

      field_grad[d][d] += qj*Br;
      field_grad[d][0] += (-qj*Cr*rij[d]*rij[0]);
      field_grad[d][1] += (-qj*Cr*rij[d]*rij[1]);
      field_grad[d][2] += (-qj*Cr*rij[d]*rij[2]);
    }
  }

  { //dipole interactions
    //charge
    energy += (qi*muj_r - qj*mui_r)*Br;
    for(int d = 0; d < DIM; d++){
      force[d] += (Cr*(qi*muj_r - qj*mui_r)*rij[d]
                   - Br*(qi*muj[d] - qj*mui[d]));
    }

    //dipole
    energy += (Br*mui_muj - Cr*mui_r*muj_r);
    for(int d = 0; d < DIM; d++){
      force[d]  += (Cr*(mui_muj*rij[d] + muj_r*mui[d] + mui_r*muj[d])
                    - Dr*mui_r*muj_r*rij[d]);
      field[d] += (-Br*muj[d] + Cr*rij[d]*muj_r);
      
      field_grad[d][d] += (Cr*muj_r);
      field_grad[d][0] += (Cr*rij[d]*muj[0] + Cr*muj[d]*rij[0] - Dr*muj_r*rij[d]*rij[0]);
      field_grad[d][1] += (Cr*rij[d]*muj[1] + Cr*muj[d]*rij[1] - Dr*muj_r*rij[d]*rij[1]);
      field_grad[d][2] += (Cr*rij[d]*muj[2] + Cr*muj[d]*rij[2] - Dr*muj_r*rij[d]*rij[2]);
    }
  }

  { //torque on dipoles
    torque[0] = mui[1]*field[2] - mui[2]*field[1];
    torque[1] = mui[2]*field[0] - mui[0]*field[2];
    torque[2] = mui[0]*field[1] - mui[1]*field[0];
  }
}
void ewald_minimum_image(double & energy,
                         double **force,
                         double **torque,
                         double **field,
			 double ***field_grad,
                         const int &Nparticles,
                         const parallelepiped &rcell,
			 int const* gid,
                         double const* const* r,
                         double const* q,
                         double const* const* mu,
                         FILE* fout,
                         const char* save_buffer
                         ){
  double dmy_energy;
  double dmy_force[DIM], dmy_torque[DIM], dmy_field[DIM], dmy_field_grad[DIM][DIM];
  double rij_mi[DIM];
  for(int i = 0; i < Nparticles; i++){
    const int    iid  = gid[i];
    const double qi   = q[i];
    const double* mui = mu[i];

    for(int j = 0; j < Nparticles; j++){
      const int   jid   = gid[j];
      const double qj   = q[j];
      const double* muj = mu[j];

      if(iid != jid){
        rcell.distance_MI(r[i], r[j], rij_mi);
        pair_interaction(rij_mi, qi, qj, mui, muj,
                         dmy_energy, dmy_force, dmy_torque, dmy_field, dmy_field_grad);
        
        energy += (dmy_energy/2.0);
        for(int d = 0; d < DIM; d++){
          force[i][d]  += dmy_force[d];
          torque[i][d] += dmy_torque[d];
          field[i][d]  += dmy_field[d];

          field_grad[i][d][0] += dmy_field_grad[d][0];
          field_grad[i][d][1] += dmy_field_grad[d][1];
          field_grad[i][d][2] += dmy_field_grad[d][2];
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
      fprintf(fsave, "%20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E\n",
              force[i][0], force[i][1], force[i][2],
              torque[i][0], torque[i][1], torque[i][2],
              field[i][0], field[i][1], field[i][2],
              field_grad[i][0][0], field_grad[i][0][1], field_grad[i][0][2],
              field_grad[i][1][0], field_grad[i][1][1], field_grad[i][1][2],
              field_grad[i][2][0], field_grad[i][2][1], field_grad[i][2][2]
              );
    }
    fclose(fsave);
  }
}
void ewald_direct_sum(double &energy, 
                      double **force, 
                      double **torque, 
                      double **field, 
		      double ***field_grad,
                      const int &ncut,
                      const int &Nparticles,
                      const parallelepiped &rcell,
		      int const* gid,
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
  double*** shell_field_grad = (double***)alloc_3d_double(Nparticles, DIM, DIM);
  double shell_energy;
  double dmy_force[DIM];
  double dmy_torque[DIM];
  double dmy_field[DIM];
  double dmy_field_grad[DIM][DIM];
  double lbox[DIM];
  double dmy_energy;
  double rij[DIM];
  clock_t start_t, end_t;
  double cpu_t;


  // generate spherical shells of periodic cells 
  rcell.get_lengths(lbox[0], lbox[1], lbox[2]);
  nshell_list(ncut, lbox[0], lbox[1], lbox[2], 
	      lmax, mmax, nmax, ntot, nmag, nlist, fout);

  start_t = clock();
  energy = shell_energy = 0.0;
  for(int i = 0; i < Nparticles; i++){
    for(int d = 0; d < DIM; d++){
      force[i][d] = shell_force[i][d] = 0.0;
      torque[i][d] = shell_torque[i][d] = 0.0;
      field[i][d] = shell_field[i][d] = 0.0;

      field_grad[i][d][0] = shell_field_grad[i][d][0] = 0.0;
      field_grad[i][d][1] = shell_field_grad[i][d][1] = 0.0;
      field_grad[i][d][2] = shell_field_grad[i][d][2] = 0.0;
    }
  }

  for(int nshell = 0; nshell < (int)nlist.size(); nshell++){
    vector< vector<int> > &shell = nlist.at(nshell);

    for(int i = 0; i < Nparticles; i++){
      shell_energy = 0.0;
      for(int d = 0; d < DIM; d++){
	shell_force[i][d] = 0.0;
	shell_torque[i][d] = 0.0;
	shell_field[i][d] = 0.0;
	
        shell_field_grad[i][d][0] = 0.0;
        shell_field_grad[i][d][1] = 0.0;
        shell_field_grad[i][d][2] = 0.0;
      }
    }

    for(int ncell = 0; ncell < shell.size(); ncell++){
      vector<int> &cell = shell.at(ncell);
      icell[0] = cell.at(0);
      icell[1] = cell.at(1);
      icell[2] = cell.at(2);

#pragma omp parallel for schedule(dynamic, 1) private(rij, dmy_energy, dmy_force, dmy_torque, dmy_field, dmy_field_grad) \
  reduction(+:shell_energy)
      for(int i = 0; i < Nparticles; i++){
	const int    iid  = gid[i];
        const double qi   = q[i];
        const double* ri  = r[i];
	const double* mui = mu[i];
	
	for(int j = 0; j < Nparticles; j++){
	  const int    jid  = gid[j];
          const double qj   = q[j];
          const double* rj  = r[j];
	  const double* muj = mu[j];

	  if(!(iid == jid && ll == 0 && mm == 0 && nn == 0)){
            
	    for(int d = 0; d < DIM; d++){
	      rij[d] = (ri[d] - rj[d]) + (double)icell[d]*lbox[d];
	    }
            pair_interaction(rij, qi, qj, mui, muj, 
                             dmy_energy, dmy_force, dmy_torque, dmy_field, dmy_field_grad);
	    
	    shell_energy += (dmy_energy/2.0);
	    for(int d = 0; d < DIM; d++){
#pragma omp atomic
	      shell_force[i][d] += dmy_force[d];
#pragma omp atomic
	      shell_torque[i][d] += dmy_torque[d];
#pragma omp atomic
	      shell_field[i][d] += dmy_field[d];
#pragma omp atomic
              shell_field_grad[i][d][0] += dmy_field_grad[d][0];
#pragma omp atomic
              shell_field_grad[i][d][1] += dmy_field_grad[d][1];
#pragma omp atomic
              shell_field_grad[i][d][2] += dmy_field_grad[d][2];
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

        field_grad[i][d][0] += shell_field_grad[i][d][0];
        field_grad[i][d][1] += shell_field_grad[i][d][1];
        field_grad[i][d][2] += shell_field_grad[i][d][2];
      }
      if(i == 0 && print_gold){
	fprintf(fout, 
		"%3d %10d %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
		nshell, (int) shell.size(),
		shell_energy, energy,
		force[0][0], force[0][1], force[0][2],
		torque[0][0], torque[0][1], torque[0][2],
		field[0][0], field[0][1], field[0][2],
		field_grad[0][0][0], field_grad[0][1][1], 
                field_grad[0][0][0] + field_grad[0][1][1] + field_grad[0][2][2],
                SQ(field_grad[0][0][1] - field_grad[0][1][0]) +
                SQ(field_grad[0][0][2] - field_grad[0][2][0]) +
                SQ(field_grad[0][1][2] - field_grad[0][2][1])
                );
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
      fprintf(fsave, "%20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E \n", 
              force[i][0], force[i][1], force[i][2],
              torque[i][0], torque[i][1], torque[i][2],
              field[i][0], field[i][1], field[i][2],
	      field_grad[i][0][0], field_grad[i][0][1], field_grad[i][0][2],
	      field_grad[i][1][0], field_grad[i][1][1], field_grad[i][1][2],
	      field_grad[i][2][0], field_grad[i][2][1], field_grad[i][2][2]
              );
    }
    fclose(fsave);
  }
  end_t = clock();
  cpu_t = ((double)end_t - start_t)/CLOCKS_PER_SEC;
  fprintf(fout, "\tExecution Time: %12.5f \n", cpu_t);

  free_2d_double(shell_force);
  free_2d_double(shell_torque);
  free_2d_double(shell_field);
  free_3d_double(shell_field_grad);
}

void ewald_direct_sum_naive(double &energy, 
                            double **force, 
                            double **torque, 
                            double **field,
			    double ***field_grad,
                            const int &lmax, 
                            const int &mmax, 
                            const int &nmax,
                            const int &Nparticles,
                            const parallelepiped &rcell,
			    int const* gid,
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
  double dmy_field_grad[DIM][DIM];
  double dmy_energy;
  double rij[DIM];

  energy = 0.0;
  for(int i = 0; i < Nparticles; i++){
    for(int d = 0; d < DIM; d++){
      force[i][d] = 0.0;
      torque[i][d] = 0.0;
      field[i][d] = 0.0;

      field_grad[i][d][0] = 0.0;
      field_grad[i][d][1] = 0.0;
      field_grad[i][d][2] = 0.0;
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
	  const int iid     = gid[i];
          const double qi   = q[i];
          const double* ri  = r[i];
	  const double* mui = mu[i];

	  for(int j = 0; j < Nparticles; j++){
	    const int jid     = gid[j];
            const double qj   = q[j];
            const double* rj  = r[j];
	    const double* muj = mu[j];
            
	    if(!(iid == jid && ll == 0 && mm == 0 && nn == 0)){
	      for(int d = 0; d < DIM; d++){
		rij[d] = (ri[d] - rj[d]) + (double)icell[d]*lbox[d];
	      }
              pair_interaction(rij, qi, qj, mui, muj,
                               dmy_energy, dmy_force, dmy_torque,dmy_field, dmy_field_grad);

	      energy+= (dmy_energy/2.0);
	      for(int d = 0; d < DIM; d++){
		force[i][d] += dmy_force[d];
		torque[i][d] += dmy_torque[d];
		field[i][d] += dmy_field[d];

                field_grad[i][d][0] += dmy_field_grad[d][0];
                field_grad[i][d][1] += dmy_field_grad[d][1];
                field_grad[i][d][2] += dmy_field_grad[d][2];
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
      fprintf(fsave, "%20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E\n", 
              force[i][0], force[i][1], force[i][2],
              torque[i][0], torque[i][1], torque[i][2],
              field[i][0], field[i][1], field[i][2],
              field_grad[i][0][0], field_grad[i][0][1], field_grad[i][0][2],
              field_grad[i][1][0], field_grad[i][1][1], field_grad[i][1][2],
              field_grad[i][2][0], field_grad[i][2][1], field_grad[i][2][2]
              );
    }
    fclose(fsave);
  }
}
