#include "ewald_gold.h"

// Compute pair interaction energy, as well as
// force, torque, and efield on particle i due to particle j
inline void pair_interaction(const double rij[DIM],
                             const double &qi, const double &qj,
                             const double mui[DIM], const double muj[DIM],
                             const double thetai[DIM*DIM], const double thetaj[DIM*DIM],
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

  double tr_thetai = thetai[0] + thetai[4] + thetai[8];
  double tr_thetaj = thetaj[0] + thetaj[4] + thetaj[8];
  double r_thetai[DIM], r_thetaj[DIM];
  double thetai_r[DIM], thetaj_r[DIM];
  double sym_thetai_r[DIM], sym_thetaj_r[DIM];
  v_M_prod(r_thetai, rij, thetai);
  v_M_prod(r_thetaj, rij, thetaj);
  M_v_prod(thetai_r, thetai, rij);
  M_v_prod(thetaj_r, thetaj, rij);
  for(int d = 0; d < DIM; d++){
    sym_thetai_r[d] = r_thetai[d] + thetai_r[d];
    sym_thetaj_r[d] = r_thetaj[d] + thetaj_r[d];
  }

  double tr_thetai_thetaj = 0.0;
  double tr_thetai_ttthetaj = 0.0;
  for(int i = 0; i < DIM; i++){
    tr_thetai_thetaj += (thetai[i*DIM]*thetaj[i] + thetai[i*DIM+1]*thetaj[DIM+i] + thetai[i*DIM+2]*thetaj[2*DIM+i]);
    tr_thetai_ttthetaj += (thetai[i*DIM]*thetaj[i*DIM] + thetai[i*DIM+1]*thetaj[i*DIM+1] + thetai[i*DIM+2]*thetaj[i*DIM+2]);
  }

  double r_thetai_r = rij[0]*thetai_r[0] + rij[1]*thetai_r[1] + rij[2]*thetai_r[2];
  double r_thetaj_r = rij[0]*thetaj_r[0] + rij[1]*thetaj_r[1] + rij[2]*thetaj_r[2];


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

  { //quadrupole interactions

    //charge
    energy += (-Br*(qi*tr_thetaj + qj*tr_thetai) + Cr*(qi*r_thetaj_r + qj*r_thetai_r) ) / 3.0;
    for(int d = 0; d < DIM; d++){
      force[d] += (Dr*rij[d]*(qi*r_thetaj_r + qj*r_thetai_r)
                   - Cr*(rij[d]*(qi*tr_thetaj + qj*tr_thetai) + (qi*sym_thetaj_r[d] + qj*sym_thetai_r[d]) ) )/3.0;
    }
    
    //dipole
    energy += (Cr*((mui_r*tr_thetaj - muj_r*tr_thetai)
                   +(mui[0]*sym_thetaj_r[0] + mui[1]*sym_thetaj_r[1] + mui[2]*sym_thetaj_r[2])
                   -(muj[0]*sym_thetai_r[0] + muj[1]*sym_thetai_r[1] + muj[2]*sym_thetai_r[2]))
               -Dr*(mui_r*r_thetaj_r - muj_r*r_thetai_r))/3.0;
    for(int d = 0; d < DIM; d++){
      force[d] += (-Cr*((mui[d]*tr_thetaj - muj[d]*tr_thetai)
                        + ((thetaj[d*DIM] + thetaj[d])*mui[0]
                           + (thetaj[d*DIM+1] + thetaj[DIM+d])*mui[1]
                           + (thetaj[d*DIM+2] + thetaj[2*DIM+d])*mui[2])
                        - ((thetai[d*DIM] + thetai[d])*muj[0]
                           + (thetai[d*DIM+1] + thetai[DIM+d])*muj[1]
                           + (thetai[d*DIM+2] + thetai[2*DIM+d])*muj[2])
                        )
                   +Dr*(rij[d]*((mui_r*tr_thetaj - muj_r*tr_thetai)
                               +(mui[0]*sym_thetaj_r[0] + mui[1]*sym_thetaj_r[1] + mui[2]*sym_thetaj_r[2])
                               -(muj[0]*sym_thetai_r[0] + muj[1]*sym_thetai_r[1] + muj[2]*sym_thetai_r[2]))
                       + (mui_r*sym_thetaj_r[d] - muj_r*sym_thetai_r[d])
                       + (r_thetaj_r*mui[d] - r_thetai_r*muj[d])
                       )
                   -Er*(rij[d]*(mui_r*r_thetaj_r - muj_r*r_thetai_r)) )/3.0;
    }

    //quadrupole
    energy += (Cr*(tr_thetai*tr_thetaj + tr_thetai_thetaj + tr_thetai_ttthetaj)
               -Dr*((r_thetai_r*tr_thetaj + r_thetaj_r*tr_thetai)
                    +(sym_thetai_r[0]*sym_thetaj_r[0] + sym_thetai_r[1]*sym_thetaj_r[1] + sym_thetai_r[2]*sym_thetaj_r[2]))
               +Er*r_thetai_r*r_thetaj_r
               )/9.0;
    for(int d = 0; d < DIM; d++){
      force[d] += (Dr*(rij[d]*(tr_thetai*tr_thetaj + tr_thetai_thetaj + tr_thetai_ttthetaj)
                        + (tr_thetai*sym_thetaj_r[d] + tr_thetaj*sym_thetai_r[d])
                        + ((thetai[d*DIM] + thetai[d])*sym_thetaj_r[0] 
                           + (thetai[d*DIM+1] + thetai[DIM+d])*sym_thetaj_r[1] 
                           + (thetai[d*DIM+2] + thetai[2*DIM+d])*sym_thetaj_r[2])
                        + ((thetaj[d*DIM] + thetaj[d])*sym_thetai_r[0]
                           + (thetaj[d*DIM+1] + thetaj[DIM+d])*sym_thetai_r[1]
                           + (thetaj[d*DIM+2] + thetaj[2*DIM+d])*sym_thetai_r[2]
                           )
                        )
                   -Er*(rij[d]*(tr_thetai*r_thetaj_r + tr_thetaj*r_thetai_r
                                +(sym_thetai_r[0]*sym_thetaj_r[0] + sym_thetai_r[1]*sym_thetaj_r[1] + sym_thetai_r[2]*sym_thetaj_r[2]))
                        + sym_thetai_r[d]*r_thetaj_r
                        + sym_thetaj_r[d]*r_thetai_r)
                   +Fr*rij[d]*r_thetai_r*r_thetaj_r
                   )/9.0;

      field[d] += (Dr*rij[d]*r_thetaj_r - Cr*(rij[d]*tr_thetaj + sym_thetaj_r[d]))/3.0;
      
      field_grad[d][d] += (-Cr*tr_thetaj + Dr*r_thetaj_r)/3.0;
      field_grad[d][0] += (-Cr*(thetaj[d*DIM] + thetaj[d])
                           +Dr*(rij[d]*rij[0]*tr_thetaj + rij[d]*sym_thetaj_r[0] + sym_thetaj_r[d]*rij[0])
                           -Er*(rij[d]*rij[0]*r_thetaj_r)
                           )/3.0;
      field_grad[d][1] += (-Cr*(thetaj[d*DIM+1] + thetaj[DIM+d])
                           +Dr*(rij[d]*rij[1]*tr_thetaj + rij[d]*sym_thetaj_r[1] + sym_thetaj_r[d]*rij[1])
                           -Er*(rij[d]*rij[1]*r_thetaj_r)
                           )/3.0;
      field_grad[d][2] += (-Cr*(thetaj[d*DIM+2] + thetaj[2*DIM+d])
                           +Dr*(rij[d]*rij[2]*tr_thetaj + rij[d]*sym_thetaj_r[2] + sym_thetaj_r[d]*rij[2])
                           -Er*(rij[d]*rij[2]*r_thetaj_r)
                           )/3.0;
    }
  }

  { //torque on dipoles
    torque[0] = mui[1]*field[2] - mui[2]*field[1];
    torque[1] = mui[2]*field[0] - mui[0]*field[2];
    torque[2] = mui[0]*field[1] - mui[1]*field[0];
  }
  { //torque on quadrupoles (Jackson, p. 171)
    // \tau_i = 1/3 \epsilon_{ijk} Q_{jl} E_{kl}
    torque[0] += ((thetai[DIM]*field_grad[0][2] + thetai[DIM+1]*field_grad[1][2] + thetai[DIM+2]*field_grad[2][2])
                  -(thetai[2*DIM]*field_grad[0][1] + thetai[2*DIM+1]*field_grad[1][1] + thetai[2*DIM+2]*field_grad[2][1]))/3.0;
    
    torque[1] += ((thetai[2*DIM]*field_grad[0][0] + thetai[2*DIM+1]*field_grad[1][0] + thetai[2*DIM+2]*field_grad[2][0])
                  -(thetai[0]*field_grad[0][2] + thetai[1]*field_grad[1][2] + thetai[2]*field_grad[2][2]))/3.0;
    
    torque[2] += ((thetai[0]*field_grad[0][1] + thetai[1]*field_grad[1][1] + thetai[2]*field_grad[2][1]) 
                  -(thetai[DIM]*field_grad[0][0] + thetai[DIM+1]*field_grad[1][0] + thetai[DIM+2]*field_grad[2][0]))/3.0;
  }
}
void ewald_minimum_image(double & energy,
                         double **force,
                         double **torque,
                         double **field,
			 double ***field_grad,
                         const int &Nparticles,
                         const parallelepiped &rcell,
                         double const* const* r,
                         double const* q,
                         double const* const* mu,
                         double const* const* const* theta,
                         FILE* fout,
                         const char* save_buffer
                         ){
  double dmy_energy;
  double dmy_force[DIM], dmy_torque[DIM], dmy_field[DIM], dmy_field_grad[DIM][DIM];
  double rij_mi[DIM];
  for(int i = 0; i < Nparticles; i++){
    const double qi   = q[i];
    const double* mui = mu[i];
    const double* thetai = theta[i][0];

    for(int j = 0; j < Nparticles; j++){
      const double qj   = q[j];
      const double* muj = mu[j];
      const double* thetaj = theta[j][0];

      if(i != j){
        rcell.distance_MI(r[i], r[j], rij_mi);
        pair_interaction(rij_mi, qi, qj, mui, muj, thetai, thetaj,
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
                      double const* const* r,
                      double const* q,
                      double const* const* mu,
                      double const* const* const* theta,
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

    for(int ncell = 0; ncell < shell.size(); ncell++){
      vector<int> &cell = shell.at(ncell);
      icell[0] = cell.at(0);
      icell[1] = cell.at(1);
      icell[2] = cell.at(2);
      
      for(int i = 0; i < Nparticles; i++){
        const double qi = q[i];
        const double* ri = r[i];
	const double* mui = mu[i];
        const double* thetai = theta[i][0];
	
	for(int j = 0; j < Nparticles; j++){
          const double qj = q[j];
          const double* rj = r[j];
	  const double* muj = mu[j];
          const double* thetaj = theta[j][0];

	  if(!(i == j && ll == 0 && mm == 0 && nn == 0)){
            
	    for(int d = 0; d < DIM; d++){
	      rij[d] = (ri[d] - rj[d]) + (double)icell[d]*lbox[d];
	    }
            pair_interaction(rij, qi, qj, mui, muj, thetai, thetaj, 
                             dmy_energy, dmy_force, dmy_torque, dmy_field, dmy_field_grad);
	    
	    shell_energy += (dmy_energy/2.0);
	    for(int d = 0; d < DIM; d++){
	      shell_force[i][d] += dmy_force[d];
	      shell_torque[i][d] += dmy_torque[d];
	      shell_field[i][d] += dmy_field[d];

              shell_field_grad[i][d][0] += dmy_field_grad[d][0];
              shell_field_grad[i][d][1] += dmy_field_grad[d][1];
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
      if(i == 0){
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
      shell_energy = 0.0;
      for(int d = 0; d < DIM; d++){
	shell_force[i][d] = 0.0;
	shell_torque[i][d] = 0.0;
	shell_field[i][d] = 0.0;

        shell_field_grad[i][d][0] = 0.0;
        shell_field_grad[i][d][1] = 0.0;
        shell_field_grad[i][d][2] = 0.0;
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
                            double const* const* r,
                            double const* q,
                            double const* const* mu,
                            double const* const* const* theta,
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
          const double qi   = q[i];
          const double* ri  = r[i];
	  const double* mui = mu[i];
          const double* thetai = theta[i][0];

	  for(int j = 0; j < Nparticles; j++){
            const double qj   = q[j];
            const double* rj  = r[j];
	    const double* muj = mu[j];
            const double* thetaj = theta[j][0];
            
	    if(!(i == j && ll == 0 && mm == 0 && nn == 0)){
	      for(int d = 0; d < DIM; d++){
		rij[d] = (ri[d] - rj[d]) + (double)icell[d]*lbox[d];
	      }
              pair_interaction(rij, qi, qj, mui, muj, thetai, thetaj,
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
