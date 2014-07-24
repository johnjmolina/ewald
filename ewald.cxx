#include "ewald.h"

parallelepiped::parallelepiped(const double a[DIM], const double b[DIM], 
                               const double c[DIM]){
  /*
    Let (') denote quantities referring to new cell frame, otherwise
    assume cartesian lab frame (basis vectors e_i=e^i, such that x_i=x^i)

    E_i' = (tLambda)_{i'j}  e_j   -> basis vector
    e_i  = (tiLambda)_{ij'} E_j'    

    E^i' = (iLambda)_{i'j}  e_j   -> dual basis vector
    e_i  = (Lambda)_{ij'}   E^j'  

    x^i' = (iLambda)_{i'j}  x^j   -> contravariant vector
    x^i  = (Lambda)_{ij'}   x^j'  

    x_i' = (tLambda)_{i'j}  x^j   -> covariant vector
    x^i  = (tiLambda)_{ij'} x_j'  

    with Lambda = [E_1', E_2', E_3'] = [a, b, c]
    a,b,c are the cell edge vectors spanning the parallelepiped

    Note: valid k-vectors (in dual cell space) are of the form
          k_i' = 2*pi*n_i'     (n_i' are integers)
          k    = k_i' * E^i'
          
          in lab (cartesian coordinates) k vectors are thus
          k^i  = k_i  = (tiLambda)_{ij'} k_j'
                      = 2*pi (tiLambda)_{ij'} n_j'
          k    = k^i * e_i

          
   */
  tLambda[0][0] = a[0];
  tLambda[0][1] = a[1];
  tLambda[0][2] = a[2];
  tLambda[1][0] = b[0];
  tLambda[1][1] = b[1];
  tLambda[1][2] = b[2];
  tLambda[2][0] = c[0];
  tLambda[2][1] = c[1];
  tLambda[2][2] = c[2];
  M_trans(Lambda, tLambda);
  M_inv(iLambda, Lambda);
  M_inv(tiLambda, tLambda);
  
  //basis vectors E_i'
  e0 = tLambda[0];
  e1 = tLambda[1];
  e2 = tLambda[2];
  
  //dual basis vectors E^i'
  E0 = iLambda[0];
  E1 = iLambda[1];
  E2 = iLambda[2];
  
  //covariant metric tensor g_{i'j'}
  gg[0][0] = v_inner_prod(e0, e0);
  gg[1][1] = v_inner_prod(e1, e1);
  gg[2][2] = v_inner_prod(e2, e2);
  gg[0][1] = gg[1][0] = v_inner_prod(e0, e1);
  gg[0][2] = gg[2][0] = v_inner_prod(e0, e2);
  gg[1][2] = gg[2][1] = v_inner_prod(e1, e2);
  
  //contravariant metric tensor g^{i'j'}
  GG[0][0] = v_inner_prod(E0, E0);
  GG[1][1] = v_inner_prod(E1, E1);
  GG[2][2] = v_inner_prod(E2, E2);
  GG[0][1] = GG[1][0] = v_inner_prod(E0, E1);
  GG[0][2] = GG[2][0] = v_inner_prod(E0, E2);
  GG[1][2] = GG[2][1] = v_inner_prod(E1, E2);   
  
  //volume
  Vol = ABS(M_det(Lambda));
  iVol= 1.0/Vol;
  PI4iVol = PI4*iVol;
  PI8iVol = 2.0*PI4iVol;
  
  //cell edges
  l0 = v_norm(e0);
  l1 = v_norm(e1);
  l2 = v_norm(e2);
  
  //angle cosines
  alpha = v_inner_prod(e0, e1)/(l0 * l1);
  beta  = v_inner_prod(e0, e2)/(l0 * l2);
  gamma = v_inner_prod(e1, e2)/(l1 * l2);
  
  //perpendicular widths
  {
    double e12[DIM], e20[DIM], e01[DIM];
    v_cross(e12, e1, e2);
    v_cross(e20, e2, e0);
    v_cross(e01, e0, e1);
    
    w0 = ABS(v_inner_prod(e0, e12)) / v_norm(e12);
    w1 = ABS(v_inner_prod(e1, e20)) / v_norm(e20);
    w2 = ABS(v_inner_prod(e2, e01)) / v_norm(e01);
    wmin = MIN(w0, MIN(w1, w2));
    wmax = MAX(w0, MAX(w1, w2));
  }
}
parallelepiped::~parallelepiped(){
}

double parallelepiped::norm_co(const double co_v[DIM])const {
  return (co_v[0] * (GG[0][0]*co_v[0] + GG[0][1]*co_v[1] + GG[0][2]*co_v[2])
          + co_v[1] * (GG[1][0]*co_v[0] + GG[1][1]*co_v[1] + GG[1][2]*co_v[2])
          + co_v[2] * (GG[2][0]*co_v[0] + GG[2][1]*co_v[2] + GG[2][2]*co_v[2]));
}
double parallelepiped::norm_contra(const double contra_v[DIM])const {
  return (contra_v[0] * (gg[0][0]*contra_v[0] + gg[0][1]*contra_v[1] + gg[0][2]*contra_v[2])
          + contra_v[1] * (gg[1][0]*contra_v[0] + gg[1][1]*contra_v[1] + gg[1][2]*contra_v[2])
          + contra_v[2] * (gg[2][0]*contra_v[0] + gg[2][1]*contra_v[1] + gg[2][2]*contra_v[2]));
}
void parallelepiped::distance_MI(const double r1[DIM], const double r2[DIM], 
                                 double r12[DIM]) const{
  double s1[DIM], s2[DIM], s12[DIM];
  M_v_prod(s1, iLambda, r1);
  M_v_prod(s2, iLambda, r2);
  for(int d = 0; d < DIM; d++){
    s12[d] = s1[d] - s2[d];
    s12[d] = s12[d] - static_cast<double>(Nint(s12[d]));
  }
  M_v_prod(r12, Lambda, s12);
}

const double ewald::mu_zero[DIM] = {0.0, 0.0, 0.0};
const double ewald::theta_zero[DIM*DIM] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void ewald::copy_cell(const int &domain, int** old_cell, int** &new_cell){
  new_cell = alloc_2d_int(domain, DIM);
  for(int n = 0; n < domain; n++){
    for(int d = 0; d < DIM; d++){
      new_cell[n][d] = old_cell[n][d];
    }
  }
  free_2d_int(old_cell);
}
void ewald::copy_cell(const int &domain, double** old_cell, double** &new_cell){
  new_cell = alloc_2d_double(domain, DIM);
  for(int n = 0; n < domain; n++){
    for(int d = 0; d < DIM; d++){
      new_cell[n][d] = old_cell[n][d];
    }
  }
  free_2d_double(old_cell);
}


void ewald::init_domain_k(const double &ewald_delta, const double &ewald_conv){
  double* tiH[DIM] = {cell->tiLambda[0], cell->tiLambda[1], 
                      cell->tiLambda[2]};
  double* GG[DIM] = {cell->GG[0], cell->GG[1], cell->GG[2]};

  //set preliminary size from convergence parameters
  k2max = 4.0*eta2*(-log(ewald_delta));
  int nkmax = 0;
  for(int d = 0; d < DIM; d++){
    nkmax = MAX(nkmax, (int)( sqrt(k2max/GG[d][d]) / PI2) ) ;
  }
  fprintf(stderr, "# Ewald nkmax : %d %g %g\n", nkmax, k2max, GG[0][0]);
  nkmax += 1;
  k2max = k2max*ewald_conv;
  
  int n_mesh = 4*POW3((nkmax + 1));
  int **dmy_cell = alloc_2d_int(n_mesh, DIM);
  double **dmy_k = alloc_2d_double(n_mesh, DIM);
  
  //k-vectors in oblique coordinates
  double k_co[DIM];
  
  //k-vectors in cartesian coordinates
  double k0_0, k0_1, k0_2;
  double k1_0, k1_1, k1_2;
  double k2_0, k2_1, k2_2;
  
  double kksq;
  int mesh, max_l, max_m, max_n, min_l, min_m, min_n;
  mesh = 0;
  max_l = max_m = max_n = 0;
  min_l = min_m = 0;
  min_n = 1;

  //Precompute all k-vectors with k2 <= k2max
  //Only consider kx>0 half of plane
  for(int ll = min_l; ll <= nkmax; ll++){
    k_co[0] = PI2*static_cast<double>(ll);
    k0_0 = k_co[0]*tiH[0][0];
    k1_0 = k_co[0]*tiH[1][0];
    k2_0 = k_co[0]*tiH[2][0];
    
    for(int mm = min_m; mm <= nkmax; mm++){
      k_co[1] = PI2*static_cast<double>(mm);
      k0_1 = k0_0 + k_co[1]*tiH[0][1];
      k1_1 = k1_0 + k_co[1]*tiH[1][1];
      k2_1 = k2_0 + k_co[1]*tiH[2][1];
      
      for(int nn = min_n; nn <= nkmax; nn++){
        k_co[2] = PI2*static_cast<double>(nn);

        //cartesian components of k vectors
        k0_2 = k0_1 + k_co[2]*tiH[0][2];
        k1_2 = k1_1 + k_co[2]*tiH[1][2];
        k2_2 = k2_1 + k_co[2]*tiH[2][2];
        
        kksq = SQ(k0_2) + SQ(k1_2) + SQ(k2_2);
        if(kksq <= k2max){
          assert(mesh < n_mesh);
          if(!equal_tol(kksq, cell->norm_co(k_co), LARGE_TOL_MP)){
            fprintf(stderr, "$ %d %d %d %.16f %.16f\n",
                    ll, mm, nn, kksq, cell->norm_co(k_co)
                    );
            assert(false);
          }
          
          //covariant components in cell frame
          dmy_cell[mesh][0] = ll;
          dmy_cell[mesh][1] = mm;
          dmy_cell[mesh][2] = nn;
          
          //contravariant components in lab frame
          dmy_k[mesh][0] = k0_2;
          dmy_k[mesh][1] = k1_2;
          dmy_k[mesh][2] = k2_2;
          
          max_l = MAX(max_l, ABS(ll));
          max_m = MAX(max_m, ABS(mm));
          max_n = MAX(max_n, ABS(nn));
          mesh++;
        }
      }//nn
      min_n = -nkmax;
    }//mm
    min_m = -nkmax;
  }//ll
  
  ewald_domain = mesh;
  copy_cell(ewald_domain, dmy_cell, ewald_cell);
  copy_cell(ewald_domain, dmy_k, ewald_k);
  kmax_l = max_l;
  kmax_m = max_m;
  kmax_n = max_n;
  
  coskr_l = alloc_2d_double(kmax_l + 1, nump);
  sinkr_l = alloc_2d_double(kmax_l + 1, nump);

  coskr_m = alloc_2d_double(kmax_m + 1, nump);
  sinkr_m = alloc_2d_double(kmax_m + 1, nump);
  
  coskr_n = alloc_2d_double(kmax_n + 1, nump);
  sinkr_n = alloc_2d_double(kmax_n + 1, nump);

  coskr_lm = alloc_1d_double(nump);
  sinkr_lm = alloc_1d_double(nump);

  coskr = alloc_1d_double(nump);
  sinkr = alloc_1d_double(nump);


  fprintf(stderr, "# Number of k-points: %d\n", ewald_domain);
  fprintf(stderr, "# k cut (kx, ky, kz): %d %d %d\n", 
          kmax_l, kmax_m, kmax_n);
  {
    FILE* kout = filecheckopen("k_vec.dat", "w");
    fprintf(kout, "# %d %12.6f %12.6f\n", ewald_domain, 
            v_norm(ewald_k[0]), v_norm(ewald_k[ewald_domain-1]));
    for(int i = 0; i < ewald_domain; i++){
      fprintf(kout, "%10d %5d %5d %5d %16.6f\n",
              i, ewald_cell[i][0], ewald_cell[i][1], ewald_cell[i][2],
              v_norm(ewald_k[i])
            );
    }
    fclose(kout);
  }
}

void ewald::free_domain_k(){
  free_2d_int(ewald_cell);
  free_2d_double(ewald_k);

  free_2d_double(coskr_l);
  free_2d_double(sinkr_l);
  free_2d_double(coskr_m);
  free_2d_double(sinkr_m);
  free_2d_double(coskr_n);
  free_2d_double(sinkr_n);
  free_1d_double(coskr_lm);
  free_1d_double(sinkr_lm);
  free_1d_double(coskr);
  free_1d_double(sinkr);
}

ewald::ewald(parallelepiped *_cell,
             const double &ewald_alpha, 
             const double &ewald_epsilon, 
             const double &ewald_delta, 
             const double &ewald_conv,
             const int &num_particles,
             const bool& with_charge,
             const bool& with_dipole,
             const bool& with_quadrupole){
  cell = _cell;

  {//particle parameters
  nump = num_particles;
    group  = alloc_1d_int(nump);
    for(int i = 1; i <= nump; i++){
      group[i] = -i;
    }
  }
  
  {//ewald parameters
    rcut = (cell->wmin)/2.0;
    r2max= rcut*rcut;
    eta  = ewald_alpha/(rcut*2.0);
    fprintf(stderr, "#ewald alpha = %12.6E\n", eta);
    eta2 = SQ(eta);
    eta3 = eta2*eta;
    eta_exp = -1.0 / (4.0 * eta2);
    epsilon_bnd = ewald_epsilon;

    TINFOIL = (epsilon_bnd < 0 ? true : false);
    CHARGE = (with_charge ? true : false);
    DIPOLE = (with_dipole ? true : false);
    QUADRUPOLE = (with_quadrupole ? true : false);

    this->init_domain_k(ewald_delta, ewald_conv);
  }
}
ewald::~ewald(){
  free_1d_int(group);
  this->free_domain_k();
}

void ewald::add_group(const int &gid, const int &num_elem, const int* pid){
  for(int i = 0; i < num_elem; i++){
    group[pid[i]] = gid;
  }
}

void ewald::reset(double *force, 
		  double *torque, 
		  double *efield,
		  double *efield_grad){
#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < nump; i++){
    const int ii = i * DIM;
    const int iii= ii * DIM;

    force[ii] = force[ii+1] = force[ii+2] = 0.0;
    torque[ii] = torque[ii+1] = torque[ii+2] = 0.0;
    efield[ii] = efield[ii+1] = efield[ii+2] = 0.0;

    efield_grad[iii] = efield_grad[iii+1] = efield_grad[iii+2] = 0.0;
    efield_grad[iii+3] = efield_grad[iii+4] = efield_grad[iii+5] = 0.0;
    efield_grad[iii+6] = efield_grad[iii+7] = efield_grad[iii+8] = 0.0;
  }
}

void ewald::reset_boundary(const double &ewald_epsilon){
  epsilon_bnd = ewald_epsilon;
  TINFOIL = (epsilon_bnd < 0 ? true : false);
}

void ewald::compute_self(double &energy, 
			 double* force, 
			 double* torque, 
			 double* efield,
			 double* efield_grad,
                         double const* r, 
			 double const* q, 
			 double const* mu, 
			 double const* theta) const{

  const double eta_factor = eta*iRoot_PI;
  const double eta3_factor = 4.0/3.0*eta2*eta_factor;
  const double eta5_factor = 2.0/5.0*eta2*eta3_factor;

  if(CHARGE){
    double dmy_energy = 0.0;
    double dmy_grad   = 0.0;
#pragma omp parallel for schedule(dynamic, 1) private(dmy_grad) \
  reduction(+:dmy_energy)
    for(int i = 0; i < nump; i++){
      const int ii = i * DIM;
      const int iii= ii * DIM;
      dmy_energy += q[i]*q[i];

      dmy_grad = (-eta3_factor*q[i]);
      efield_grad[iii]   += dmy_grad;  //xx
      efield_grad[iii+4] += dmy_grad;  //yy
      efield_grad[iii+8] += dmy_grad;  //zz
    }
    energy += (-eta_factor*dmy_energy);
  }

  if(DIPOLE){
    double dmy_energy = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_energy)
    for(int i = 0; i < nump; i++){
      const int ii = i * DIM;
      const double* mui = &mu[ii];
      dmy_energy   += mui[0]*mui[0] + mui[1]*mui[1] + mui[2]*mui[2];
      efield[ii]   += eta3_factor*mui[0];
      efield[ii+1] += eta3_factor*mui[1];
      efield[ii+2] += eta3_factor*mui[2];
    }
    energy += (-eta3_factor*dmy_energy/2.0);
  }

  if(QUADRUPOLE){
    double dmy_energy_q = 0.0;
    double dmy_energy_theta = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_energy_q,dmy_energy_theta)
    for(int i = 0; i < nump; i++){
      const int iii = i*DIM*DIM;
      const double qi = (CHARGE ? q[i] : 0.0);
      const double* thetai = &theta[iii];

      dmy_energy_q += qi*(thetai[0] + thetai[4] + thetai[8]);
      for(int j = 0; j < DIM; j++){
        const int jj = j * DIM;
        dmy_energy_theta += (thetai[jj]*thetai[jj] + thetai[jj+1]*thetai[jj+1] + thetai[jj+2]*thetai[jj+2]);
      }
    }
    energy += (-eta3_factor*dmy_energy_q - eta5_factor*dmy_energy_theta)/3.0;
  }
}

void ewald::compute_surface(double &energy, 
			    double* force, 
			    double* torque, 
			    double* efield,
			    double* efield_grad,
                            double const* r, 
			    double const* q, 
			    double const* mu, 
			    double const* theta) const{
  if(!TINFOIL){
    
    double sum_qr0, sum_qr1, sum_qr2;
    sum_qr0 = sum_qr1 = sum_qr2 = 0.0;
    if(CHARGE){
#pragma omp parallel for schedule(dynamic, 1) reduction(+:sum_qr0, sum_qr1, sum_qr2)
      for(int i = 0; i < nump; i++){
        const double* ri = &r[i*DIM];
        sum_qr0 += q[i]*ri[0];
        sum_qr1 += q[i]*ri[1];
        sum_qr2 += q[i]*ri[2];
      }
    }

    double sum_mu0, sum_mu1, sum_mu2;
    sum_mu0 = sum_mu1 = sum_mu2 = 0.0;
    if(DIPOLE){
#pragma omp parallel for schedule(dynamic, 1) reduction(+:sum_mu0, sum_mu1, sum_mu2)
      for(int i = 0; i < nump; i++){
        const double* mui = &mu[i*DIM];
        sum_mu0 += mui[0];
        sum_mu1 += mui[1];
        sum_mu2 += mui[2];
      }
    }

    double dmy_factor = (2.0*M_PI / (2.0*epsilon_bnd + 1.0)) * (cell->iVol);
    double sum_qr_mu0= sum_qr0 + sum_mu0;
    double sum_qr_mu1= sum_qr1 + sum_mu1;
    double sum_qr_mu2= sum_qr2 + sum_mu2;
    energy += dmy_factor*(sum_qr_mu0*sum_qr_mu0 + sum_qr_mu1*sum_qr_mu1 + sum_qr_mu2*sum_qr_mu2);

    dmy_factor *= (-2.0);
#pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < nump; i++){
      const int ii = i * DIM;
      const int iii = ii * DIM;
      const double  qi  = (CHARGE ? q[i] : 0.0);
      const double* mui = (DIPOLE ? &mu[ii] : mu_zero);
      const double* thetai = (QUADRUPOLE ? &theta[iii] : theta_zero);


      efield[ii]   += dmy_factor*sum_qr_mu0;
      efield[ii+1] += dmy_factor*sum_qr_mu1;
      efield[ii+2] += dmy_factor*sum_qr_mu2;

      if(CHARGE){
        force[ii]   += dmy_factor*qi*sum_qr_mu0;
        force[ii+1] += dmy_factor*qi*sum_qr_mu1;
        force[ii+2] += dmy_factor*qi*sum_qr_mu2;
      }
      
      if(DIPOLE){
        torque[ii]   += dmy_factor*(mui[1]*sum_qr_mu2 - mui[2]*sum_qr_mu1);
        torque[ii+1] += dmy_factor*(mui[2]*sum_qr_mu0 - mui[0]*sum_qr_mu2);
        torque[ii+2] += dmy_factor*(mui[0]*sum_qr_mu1 - mui[1]*sum_qr_mu0);
      }
    }
  }
}

void ewald::compute_r(double &energy, 
		      double* force, 
		      double* torque, 
		      double* efield,
		      double* efield_grad,
                      double const* r, 
		      double const* q, 
		      double const* mu, 
		      double const* theta) const{
  
  double dmy_energy = 0.0;
  
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_energy)
  for(int j = 1; j < nump; j++){
    const int jj = j*DIM;
    const int jjj= jj*DIM;
    const int jid= group[j];
    
    double rij[DIM];
    double erfc_ewald;
    double drij, drij2;
    double Br, Cr, Dr, Er;
    double dmy_0, dmy_1, dmy_2, dmy_3;
    
    double dmy_force[DIM];
    double dmy_efieldi[DIM], dmy_efieldj[DIM];
    double dmy_efieldi_dx[DIM], dmy_efieldi_dy[DIM], dmy_efieldi_dz[DIM];
    double dmy_efieldj_dx[DIM], dmy_efieldj_dy[DIM], dmy_efieldj_dz[DIM];
    
    const double  qj  = (CHARGE ? q[j] : 0.0);
    const double* muj = (DIPOLE ? &mu[jj] : mu_zero);
    const double* thetaj = (QUADRUPOLE ? &theta[jjj] : theta_zero);
    for(int i = 0; i < j; i++){
      const int ii = i*DIM;
      const int iii= ii*DIM;
      const int iid= group[i];
      
      //i and j belong to different group
      if(jid != iid){
        
        cell->distance_MI(&r[ii], &r[jj], rij);
        drij = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	//i and j within cutoff distance
        if(drij < r2max){
          const double  qi  = (CHARGE ? q[i] : 0.0);
          const double* mui = (DIPOLE ? &mu[ii] : mu_zero);
          const double* thetai = (QUADRUPOLE ? &theta[iii] : theta_zero);
          dmy_force[0] = dmy_force[1] = dmy_force[2] = 0.0;
          dmy_efieldi[0] = dmy_efieldi[1] = dmy_efieldi[2] = 0.0;
          dmy_efieldj[0] = dmy_efieldj[1] = dmy_efieldj[2] = 0.0;
          
          dmy_efieldi_dx[0] = dmy_efieldi_dx[1] = dmy_efieldi_dx[2] = 0.0;
          dmy_efieldi_dy[0] = dmy_efieldi_dy[1] = dmy_efieldi_dy[2] = 0.0;
          dmy_efieldi_dz[0] = dmy_efieldi_dz[1] = dmy_efieldi_dz[2] = 0.0;
          
          dmy_efieldj_dx[0] = dmy_efieldj_dx[1] = dmy_efieldj_dx[2] = 0.0;
          dmy_efieldj_dy[0] = dmy_efieldj_dy[1] = dmy_efieldj_dy[2] = 0.0;
          dmy_efieldj_dz[0] = dmy_efieldj_dz[1] = dmy_efieldj_dz[2] = 0.0;
          
          drij = sqrt(drij);
          erfc_ewald = erfc(eta*drij);
          dmy_0 = exp(-eta2*drij*drij);
          drij = 1.0 / drij;
          drij2= drij*drij;

          dmy_0 *= 2.0*eta*iRoot_PI;
          Br = (erfc_ewald*drij + dmy_0) * drij2;
          
          dmy_0 *= 2.0*eta2;
          Cr = (3.0 * Br + dmy_0) * drij2;
          
          dmy_0 *= 2.0*eta2;
          Dr = (5.0 * Cr + dmy_0) * drij2;

          dmy_0 *= 2.0*eta2;
          Er = (7.0 * Dr + dmy_0) * drij2;
          
          if(CHARGE){
            dmy_energy += erfc_ewald*drij*qi*qj;
            
            for(int d = 0; d < DIM; d++){
              dmy_0 = Br*rij[d];
              dmy_force[d]   += qj*qi*dmy_0;
              dmy_efieldi[d] += qj*dmy_0;
              dmy_efieldj[d] += (-qi*dmy_0);
              
              dmy_0 = (-Cr*rij[d]);
              dmy_1 = qj * dmy_0;
              dmy_2 = qi * dmy_0;
              dmy_efieldi_dx[d] += (dmy_1*rij[0]);
              dmy_efieldi_dy[d] += (dmy_1*rij[1]);
              dmy_efieldi_dz[d] += (dmy_1*rij[2]);
              dmy_efieldj_dx[d] += (dmy_2*rij[0]);
              dmy_efieldj_dy[d] += (dmy_2*rij[1]);
              dmy_efieldj_dz[d] += (dmy_2*rij[2]);
            }
            
            dmy_0 = qj*Br;
            dmy_1 = qi*Br;
            dmy_efieldi_dx[0] += dmy_0;
            dmy_efieldi_dy[1] += dmy_0;
            dmy_efieldi_dz[2] += dmy_0;
            dmy_efieldj_dx[0] += dmy_1;
            dmy_efieldj_dy[1] += dmy_1;
            dmy_efieldj_dz[2] += dmy_1;
          }

          if(DIPOLE){
            double mui_r, muj_r, mui_muj;
            mui_r = v_inner_prod(mui, rij);
            muj_r = v_inner_prod(muj, rij);
            mui_muj = v_inner_prod(mui, muj);
            
            //charge - dipole
            if(CHARGE){
              dmy_energy += Br*(qi*muj_r - qj*mui_r);
              dmy_0 = Cr*(qi*muj_r - qj*mui_r);
              dmy_1 = qj*Br;
              dmy_2 = -qi*Br;
              for(int d = 0; d < DIM; d++) dmy_force[d] += (dmy_0*rij[d] + dmy_1*mui[d] + dmy_2*muj[d]);
            }
            
            //dipole - dipole          
            dmy_energy += Br*mui_muj - Cr*mui_r*muj_r;
            for(int d = 0; d < DIM; d++){
              dmy_0 = rij[d] * Cr;
              dmy_1 = muj[d] * Cr;
              dmy_2 = mui[d] * Cr;
              dmy_3 = -rij[d] * muj_r * Dr;
              
              dmy_force[d] += (dmy_0*mui_muj + (dmy_1+dmy_3)*mui_r + dmy_2*muj_r);
              dmy_efieldi[d] += (-Br*muj[d] + dmy_0*muj_r);
              dmy_efieldj[d] += (-Br*mui[d] + dmy_0*mui_r);
              
              dmy_efieldi_dx[d] += (dmy_0*muj[0] + dmy_1*rij[0] + dmy_3*rij[0]);
              dmy_efieldi_dy[d] += (dmy_0*muj[1] + dmy_1*rij[1] + dmy_3*rij[1]);
              dmy_efieldi_dz[d] += (dmy_0*muj[2] + dmy_1*rij[2] + dmy_3*rij[2]);
              
              dmy_3 = -rij[d] * mui_r * Dr;
              dmy_efieldj_dx[d] -= (dmy_0*mui[0] + dmy_2*rij[0] + dmy_3*rij[0]);
              dmy_efieldj_dy[d] -= (dmy_0*mui[1] + dmy_2*rij[1] + dmy_3*rij[1]);
              dmy_efieldj_dz[d] -= (dmy_0*mui[2] + dmy_2*rij[2] + dmy_3*rij[2]);
            }
            dmy_0 = Cr*muj_r;
            dmy_1 = -Cr*mui_r;
            dmy_efieldi_dx[0] += dmy_0;
            dmy_efieldi_dy[1] += dmy_0;
            dmy_efieldi_dz[2] += dmy_0;
            dmy_efieldj_dx[0] += dmy_1;
            dmy_efieldj_dy[1] += dmy_1;
            dmy_efieldj_dz[2] += dmy_1;
          }

          if(QUADRUPOLE){
            //quad-quad
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
            double tr_thetai_thetaj = 0;
            double tr_thetai_ttthetaj = 0;
            for(int d = 0; d < DIM; d++){
              tr_thetai_thetaj += (thetai[d*DIM]*thetaj[d] + thetai[d*DIM+1]*thetaj[DIM+d] + thetai[d*DIM+2]*thetaj[2*DIM+d]);
              tr_thetai_ttthetaj += (thetai[d*DIM]*thetaj[d*DIM] + thetai[d*DIM+1]*thetaj[d*DIM+1] + thetai[d*DIM+2]*thetaj[d*DIM+2]);
            }
            double r_thetai_r = rij[0]*thetai_r[0] + rij[1]*thetai_r[1] + rij[2]*thetai_r[2];
            double r_thetaj_r = rij[0]*thetaj_r[0] + rij[1]*thetaj_r[1] + rij[2]*thetaj_r[2];
            dmy_energy += (Cr*(tr_thetai*tr_thetaj + tr_thetai_thetaj + tr_thetai_ttthetaj)
                           -Dr*((r_thetai_r*tr_thetaj + r_thetaj_r*tr_thetai)
                                +(sym_thetai_r[0]*sym_thetaj_r[0] + sym_thetai_r[1]*sym_thetaj_r[1] + sym_thetai_r[2]*sym_thetaj_r[2]))
                           +Er*r_thetai_r*r_thetaj_r)/9.0;
          }
          
          {
            //forces
#pragma omp atomic
            force[ii]   += dmy_force[0];
#pragma omp atomic
            force[ii+1] += dmy_force[1];
#pragma omp atomic
            force[ii+2] += dmy_force[2];
            
            force[jj]   -= dmy_force[0];
            force[jj+1] -= dmy_force[1];
            force[jj+2] -= dmy_force[2];
          }
          
          {
            //electric fields
#pragma omp atomic
            efield[ii]   += dmy_efieldi[0];
#pragma omp atomic
            efield[ii+1] += dmy_efieldi[1];
#pragma omp atomic
            efield[ii+2] += dmy_efieldi[2];
            
            efield[jj]   += dmy_efieldj[0];
            efield[jj+1] += dmy_efieldj[1];
            efield[jj+2] += dmy_efieldj[2];
          }
          
          //gradient of electric fields
          {
#pragma omp atomic
            efield_grad[iii]   += dmy_efieldi_dx[0];
#pragma omp atomic
            efield_grad[iii+1] += dmy_efieldi_dx[1];
#pragma omp atomic
            efield_grad[iii+2] += dmy_efieldi_dx[2];
#pragma omp atomic
            efield_grad[iii+3] += dmy_efieldi_dy[0];
#pragma omp atomic
            efield_grad[iii+4] += dmy_efieldi_dy[1];
#pragma omp atomic
            efield_grad[iii+5] += dmy_efieldi_dy[2];
#pragma omp atomic
            efield_grad[iii+6] += dmy_efieldi_dz[0];
#pragma omp atomic
            efield_grad[iii+7] += dmy_efieldi_dz[1];
#pragma omp atomic
            efield_grad[iii+8] += dmy_efieldi_dz[2];
            
            efield_grad[jjj]   += dmy_efieldj_dx[0];
            efield_grad[jjj+1] += dmy_efieldj_dx[1];
            efield_grad[jjj+2] += dmy_efieldj_dx[2];
            efield_grad[jjj+3] += dmy_efieldj_dy[0];
            efield_grad[jjj+4] += dmy_efieldj_dy[1];
            efield_grad[jjj+5] += dmy_efieldj_dy[2];
            efield_grad[jjj+6] += dmy_efieldj_dz[0];
            efield_grad[jjj+7] += dmy_efieldj_dz[1];
            efield_grad[jjj+8] += dmy_efieldj_dz[2];
          }
          
          //torques
          if(DIPOLE){
#pragma omp atomic
            torque[ii]   += (mui[1]*dmy_efieldi[2] - mui[2]*dmy_efieldi[1]);
#pragma omp atomic
            torque[ii+1] += (mui[2]*dmy_efieldi[0] - mui[0]*dmy_efieldi[2]);
#pragma omp atomic
            torque[ii+2] += (mui[0]*dmy_efieldi[1] - mui[1]*dmy_efieldi[0]);
            
            torque[jj]   += (muj[1]*dmy_efieldj[2] - muj[2]*dmy_efieldj[1]);
            torque[jj+1] += (muj[2]*dmy_efieldj[0] - muj[0]*dmy_efieldj[2]);
            torque[jj+2] += (muj[0]*dmy_efieldj[1] - muj[1]*dmy_efieldj[0]);
          }
          if(QUADRUPOLE){
            //TODO: add torque due to gradient of electric field
          }
        }//rij < r2max
      }//jid != iid
    }//i
  }//j
  energy += dmy_energy;
}

void ewald::precompute_trig_k(double const* r){
  double* iH = (cell->iLambda)[0];
#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < nump; i++){
    const double* ri = &r[i*DIM];
    /*
      Transform UNIT k vectors to cartesian coordinates to compute k.r
      k.r = k_i . r^i = ( (tiLambda)_{ij'} k_j' ) . r^i
      = (2*pi (tiLambda)_{ij'} n_j') . r^i
      = (2*pi (tiLambda)_{ij'} delta_{j'l'}) . r^i
      = 2*pi (tiLambda)_{il'} r^i
      = 2*pi (iLambda)_{l'i} r^i
    */
    double twopibox_l = PI2*(iH[0]*ri[0] + iH[1]*ri[1] + iH[2]*ri[2]);
    double twopibox_m = PI2*(iH[3]*ri[0] + iH[4]*ri[1] + iH[5]*ri[2]);
    double twopibox_n = PI2*(iH[6]*ri[0] + iH[7]*ri[1] + iH[8]*ri[2]);
    
    
    //setup cosine / sine arrays
    coskr_l[0][i] = 1.0;
    coskr_m[0][i] = 1.0;
    coskr_n[0][i] = 1.0;

    sinkr_l[0][i] = 0.0;
    sinkr_m[0][i] = 0.0;
    sinkr_n[0][i] = 0.0;

    coskr_l[1][i] = cos(twopibox_l);
    coskr_m[1][i] = cos(twopibox_m);
    coskr_n[1][i] = cos(twopibox_n);

    sinkr_l[1][i] = sin(twopibox_l);
    sinkr_m[1][i] = sin(twopibox_m);
    sinkr_n[1][i] = sin(twopibox_n);

    for(int j = 2; j <= kmax_l; j++){
      coskr_l[j][i] = coskr_l[j-1][i] * coskr_l[1][i] - sinkr_l[j-1][i] * sinkr_l[1][i];
      sinkr_l[j][i] = sinkr_l[j-1][i] * coskr_l[1][i] + coskr_l[j-1][i] * sinkr_l[1][i];
    }

    for(int j = 2; j <= kmax_m; j++){
      coskr_m[j][i] = coskr_m[j-1][i] * coskr_m[1][i] - sinkr_m[j-1][i] * sinkr_m[1][i];
      sinkr_m[j][i] = sinkr_m[j-1][i] * coskr_m[1][i] + coskr_m[j-1][i] * sinkr_m[1][i];
    }

    for(int j = 2; j <= kmax_n; j++){
      coskr_n[j][i] = coskr_n[j-1][i] * coskr_n[1][i] - sinkr_n[j-1][i] * sinkr_n[1][i];
      sinkr_n[j][i] = sinkr_n[j-1][i] * coskr_n[1][i] + coskr_n[j-1][i] * sinkr_n[1][i];
    }
  }
}

void ewald::compute_trig_k(const int& ll, const int& mm, const int& nn){
  // Compute cos(k.r) and sin(k.r) for all particles for k = (ll, mm, nn)
  // Note: only ll > 0 half is considered
  int l = ABS(ll);
  int m = ABS(mm);
  int n = ABS(nn);
  double sign;

  // cos(p l + q m) = cos(p l) cos(q m) - sin(p l) sin(q m)
  // sin(p l + q m) = sin(p l) cos(q m) + cos(p l) sin(q m)
  sign = (mm > 0 ? 1.0 : -1.0);
#pragma omp parallel for schedule(dynamic, 1)
  for(int j = 0; j < nump; j++){
    coskr_lm[j] = coskr_l[l][j] * coskr_m[m][j] - sign * sinkr_l[l][j] * sinkr_m[m][j];
    sinkr_lm[j] = sinkr_l[l][j] * coskr_m[m][j] + sign * coskr_l[l][j] * sinkr_m[m][j];
  }

  // cos(p l + q m + r n) = cos(p l + q m) cos(r n) - sin(p l + q m) sin(r n)
  // sin(p l + q m + r n) = sin(p l + q m) cos(r n) + cos(p l + q m) sin(r n)
  sign = (nn > 0 ? 1.0 : -1.0);
#pragma omp parallel for schedule(dynamic, 1)
  for(int j = 0; j < nump; j++){
    coskr[j] = coskr_lm[j] * coskr_n[n][j] - sign * sinkr_lm[j] * sinkr_n[n][j];
    sinkr[j] = sinkr_lm[j] * coskr_n[n][j] + sign * coskr_lm[j] * sinkr_n[n][j];
  }
}

inline void ewald::compute_rho_k(double &rho_re, double &rho_im, 
                                 double const* q, double const* mu, double const* theta,
                                 double const kk[DIM], double const* coskr, double const* sinkr
                                 ) const{
  rho_re = rho_im = 0.0;
  double dmy_rho_re, dmy_rho_im;
  if(CHARGE){
    dmy_rho_re = dmy_rho_im = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_rho_re, dmy_rho_im)
    for(int j = 0; j < nump; j++){
      dmy_rho_re += q[j]*coskr[j];
      dmy_rho_im += q[j]*sinkr[j];
    }
    rho_re += dmy_rho_re;
    rho_im += dmy_rho_im;
  }
  if(DIPOLE){
    dmy_rho_re = dmy_rho_im = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_rho_re, dmy_rho_im)
    for(int j = 0; j < nump; j++){
      const double* muj = &mu[j*DIM];
      double mu_k = muj[0]*kk[0] + muj[1]*kk[1] + muj[2]*kk[2];
      dmy_rho_re += (-mu_k*sinkr[j]);
      dmy_rho_im += (mu_k*coskr[j]);
    }
    rho_re += dmy_rho_re;
    rho_im += dmy_rho_im;
  }
  if(QUADRUPOLE){
    dmy_rho_re = dmy_rho_im = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_rho_re, dmy_rho_im)
    for(int j = 0; j < nump; j++){
      const double* thetaj = &theta[j*DIM*DIM];
      double k_theta_k = -(kk[0]*(thetaj[0]*kk[0] + thetaj[1]*kk[1] + thetaj[2]*kk[2]) +
                           kk[1]*(thetaj[3]*kk[0] + thetaj[4]*kk[1] + thetaj[5]*kk[2]) +
                           kk[2]*(thetaj[6]*kk[0] + thetaj[7]*kk[1] + thetaj[8]*kk[2]))/3.0;
      dmy_rho_re += k_theta_k * coskr[j];
      dmy_rho_im += k_theta_k * sinkr[j];
    }
    rho_re += dmy_rho_re;
    rho_im += dmy_rho_im;
  }
}
void ewald::compute_k(double &energy, 
		      double* force, 
		      double* torque, 
		      double* efield, 
		      double* efield_grad,
                      double const* r, 
		      double const* q, 
		      double const* mu, 
		      double const* theta){
  double dmy_energy, dmy_force, dmy_efield, dmy_efield_grad, dmy_efield_grad0;
  dmy_energy = dmy_force = dmy_efield = dmy_efield_grad = 0.0;
  
  this->precompute_trig_k(r);
  for(int i = 0; i < ewald_domain; i++){
    //compute cos(kr) and sin(kr) terms for all particles
    this->compute_trig_k(ewald_cell[i][0], ewald_cell[i][1], ewald_cell[i][2]);

    //k vector in cartesian (lab) coordinates
    double &k0 = ewald_k[i][0];
    double &k1 = ewald_k[i][1];
    double &k2 = ewald_k[i][2];

    double rho_re, rho_im;
    this->compute_rho_k(rho_re, rho_im, q, mu, theta, ewald_k[i], coskr, sinkr);

    double kk = k0*k0 + k1*k1 + k2*k2;
    double ewald_damp  = exp(kk*eta_exp)/kk;
    dmy_energy += ewald_damp * (rho_re*rho_re + rho_im*rho_im);

    ewald_damp *= (cell->PI8iVol);

#pragma omp parallel for schedule(dynamic, 1) private(dmy_force, dmy_efield, dmy_efield_grad0, dmy_efield_grad)
    for(int j = 0; j < nump; j++){
      const int jj = j * DIM;
      const int jjj= jj * DIM;
      const double  qj  = (CHARGE ? q[j] : 0.0);
      const double* muj = (DIPOLE ? &mu[jj] : mu_zero);
      const double* thetaj = (QUADRUPOLE ? &theta[jjj] : theta_zero);
      double mu_k = (DIPOLE ? muj[0]*k0 + muj[1]*k1 + muj[2]*k2 : 0.0);
      double k_theta_k = (QUADRUPOLE ? 
                          -(k0*(thetaj[0]*k0 + thetaj[1]*k1 + thetaj[2]*k2) +
                            k1*(thetaj[3]*k0 + thetaj[4]*k1 + thetaj[5]*k2) +
                            k2*(thetaj[6]*k0 + thetaj[7]*k1 + thetaj[8]*k2))/3.0: 0.0);
      double qq_re    = -((qj + k_theta_k)*sinkr[j] + mu_k*coskr[j]);
      double qq_im    =  ((qj + k_theta_k)*coskr[j] - mu_k*sinkr[j]);

      dmy_force = -ewald_damp*(qq_re*rho_re + qq_im*rho_im);
      dmy_efield = -ewald_damp*(coskr[j]*rho_im - sinkr[j]*rho_re);
      dmy_efield_grad0 = ewald_damp*(coskr[j]*rho_re + sinkr[j]*rho_im);

      force[jj]   += dmy_force * k0;
      force[jj+1] += dmy_force * k1;
      force[jj+2] += dmy_force * k2;

      efield[jj]   += dmy_efield * k0;
      efield[jj+1] += dmy_efield * k1;
      efield[jj+2] += dmy_efield * k2;

      dmy_efield_grad = dmy_efield_grad0*k0;
      efield_grad[jjj]   += dmy_efield_grad * k0;  //xx
      efield_grad[jjj+1] += dmy_efield_grad * k1;  //xy
      efield_grad[jjj+2] += dmy_efield_grad * k2;  //xz
      
      dmy_efield_grad = dmy_efield_grad0*k1;
      efield_grad[jjj+3] += dmy_efield_grad * k0;  //yx
      efield_grad[jjj+4] += dmy_efield_grad * k1;  //yy
      efield_grad[jjj+5] += dmy_efield_grad * k2;  //yz
      
      dmy_efield_grad = dmy_efield_grad0*k2;
      efield_grad[jjj+6] += dmy_efield_grad * k0;  //zx
      efield_grad[jjj+7] += dmy_efield_grad * k1;  //zy
      efield_grad[jjj+8] += dmy_efield_grad * k2;  //zz

      if(DIPOLE){
        torque[jj]   += dmy_efield * (muj[1] * k2 - muj[2] * k1);
        torque[jj+1] += dmy_efield * (muj[2] * k0 - muj[0] * k2);
        torque[jj+2] += dmy_efield * (muj[0] * k1 - muj[1] * k0);
      }
      if(QUADRUPOLE){
        //TODO: compute efield_gradient + torque
      }
    }
  }
  energy += (cell->PI4iVol)*dmy_energy;
}

void ewald::compute(double* E_ewald, 
		    double* force, 
		    double* torque, 
		    double* efield,
		    double* efield_grad,
                    double const* r, 
		    double const* q,
		    double const* mu,
		    double const* theta,
                    const char* save_buffer
               ){
  if(!CHARGE && !DIPOLE && !QUADRUPOLE) return;
  double &e_r = E_ewald[1];
  double &e_k = E_ewald[2];
  double &e_self = E_ewald[3];
  double &e_surface = E_ewald[4];
  e_r = e_k = e_self = e_surface = 0.0;
  clock_t start_t, end_t;
  double cpu_t;
  start_t = clock();
  this->reset(force, torque, efield, efield_grad);
  this->compute_r(e_r, force, torque, efield, efield_grad, r, q, mu, theta);
  this->compute_k(e_k, force, torque, efield, efield_grad, r, q, mu, theta);
  this->compute_self(e_self, force, torque, efield, efield_grad, r, q, mu, theta);
  this->compute_surface(e_surface, force, torque, efield, efield_grad, r, q, mu, theta);
  end_t = clock();

  cpu_t = ((double)end_t - start_t)/CLOCKS_PER_SEC;
  fprintf(stderr, "\tExecution Time = %12.5f \n", cpu_t);
  E_ewald[0] = E_ewald[1] + E_ewald[2] + E_ewald[3] + E_ewald[4];

  {
    FILE* fsave = filecheckopen(save_buffer, "w");
    fprintf(fsave, "%d\n", nump);
    fprintf(fsave, "%20.12E\n", E_ewald[0]);
    for(int i = 0; i < nump; i++){
      const int ii = i*DIM;
      const int iii= ii*DIM;
      fprintf(fsave, "%20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E %20.12E\n",
              force[ii], force[ii+1], force[ii+2],
              torque[ii], torque[ii+1], torque[ii+2],
              efield[ii], efield[ii+1], efield[ii+2],
              efield_grad[iii], efield_grad[iii+1], efield_grad[iii+2],
              efield_grad[iii+3], efield_grad[iii+4], efield_grad[iii+5],
              efield_grad[iii+6], efield_grad[iii+7], efield_grad[iii+8]
              );
    }
    fclose(fsave);

    if(TINFOIL and q != NULL){
      //convert to atomic units to compare with cp2k
      // energy -> hartree
      // force  -> hartree / bohr radius
      // length -> bohr radius
      // assuming lenght units are in angstroms
      const double e_au = 0.5291772114258002;  //E_h
      const double f_au = 0.2800285208247281;  //E_h / a0
      const double r_au = 1.889726124565062;   //a0
      
      char cp2k_buffer[256];
      sprintf(cp2k_buffer, "cp2k_%s", save_buffer);
      FILE* fcp2k = filecheckopen(cp2k_buffer, "w");
      fprintf(fcp2k, "%5d\n", nump);
      fprintf(fcp2k, "i =        0, time =        0.000, E =        %12.10f\n", E_ewald[0]*e_au);
      for(int i = 0; i < nump; i++){
        const int ii = i * DIM;
        fprintf(fcp2k, "%5s   %16.10f   %16.10f   %16.10f\n",
                (q[i] > 0.0 ? "Na" : "Cl"),
                force[ii]*f_au, force[ii+1]*f_au, force[ii+2]*f_au
                );
      }
      fclose(fcp2k);
    }
  }
}
