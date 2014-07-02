#include "ewald.h"

void print_matrix(const double a0[DIM], const double a1[DIM], const double a2[DIM], const char* tag){
  fprintf(stderr, "%s : \n", tag);
  fprintf(stderr, "%.5f %.5f %.5f\n", a0[0], a0[1], a0[2]);
  fprintf(stderr, "%.5f %.5f %.5f\n", a1[0], a1[1], a1[2]);
  fprintf(stderr, "%.5f %.5f %.5f\n", a2[0], a2[1], a2[2]);
}
parallelepiped::parallelepiped(const double a[DIM], const double b[DIM], 
                               const double c[DIM]){
  //transfromation matrices
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

  
  //basis vectors
  e0 = tLambda[0];
  e1 = tLambda[1];
  e2 = tLambda[2];
  
  //dual basis vectors
  E0 = iLambda[0];
  E1 = iLambda[1];
  E2 = iLambda[2];
  
  //covariant metric tensor
  gg[0][0] = v_inner_prod(e0, e0);
  gg[1][1] = v_inner_prod(e1, e1);
  gg[2][2] = v_inner_prod(e2, e2);
  gg[0][1] = gg[1][0] = v_inner_prod(e0, e1);
  gg[0][2] = gg[2][0] = v_inner_prod(e0, e2);
  gg[1][2] = gg[2][1] = v_inner_prod(e1, e2);
  
  //contravariant metric tensor
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
  fprintf(stderr, "nmax : %d %g %g\n", nkmax, k2max, GG[0][0]);
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
          
          dmy_cell[mesh][0] = ll;
          dmy_cell[mesh][1] = mm;
          dmy_cell[mesh][2] = nn;
          
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


  fprintf(stderr, "Number of k-points: %d\n", ewald_domain);
  fprintf(stderr, "k cut (kx, ky, kz): %d %d %d\n", 
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

  //particle parameters
  nump = num_particles;
  
  {//ewald parameters
    rcut = (cell->wmin)/2.0;
    r2max= rcut*rcut;
    eta  = ewald_alpha/(rcut*2.0);
    eta2 = SQ(eta);
    eta3 = eta2*eta;
    eta_exp = -1.0 / (4.0 * eta2);
    epsilon_bnd = ewald_epsilon;

    TINFOIL = (epsilon_bnd < 0 ? true : false);
    CHARGE = false;
    DIPOLE = (with_dipole ? true : false);
    QUADRUPOLE = false;

    this->init_domain_k(ewald_delta, ewald_conv);
  }
}

void ewald::reset(double *force, double *torque, double *efield){
#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < nump; i++){
    int ii = i * DIM;
    force[ii] = force[ii+1] = force[ii+2] = 0.0;
    torque[ii] = torque[ii+1] = torque[ii+2] = 0.0;
    efield[ii] = efield[ii+1] = efield[ii+2] = 0.0;
  }
}

void ewald::reset_boundary(const double &ewald_epsilon){
  epsilon_bnd = ewald_epsilon;
  TINFOIL = (epsilon_bnd < 0 ? true : false);
}

void ewald::compute_self(double &energy, double* force, double* torque, double* efield,
                         double const* r, double const* q, double const* mu, double const* theta) const{

  double mu_factor = 4.0/3.0*eta3*SQRT_PI_inv;
  double dmy_energy = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_energy)
  for(int i = 0; i < nump; i++){
    int ii = i * DIM;
    const double* mui = &mu[ii];
    dmy_energy   += mui[0]*mui[0] + mui[1]*mui[1] + mui[2]*mui[2];
    efield[ii]   += mu_factor*mui[0];
    efield[ii+1] += mu_factor*mui[1];
    efield[ii+2] += mu_factor*mui[2];
  }
  energy += (-mu_factor*dmy_energy/2.0);
}
void ewald::compute_surface(double &energy, double* force, double* torque, double* efield,
                            double const* r, double const* q, double const* mu, double const* theta) const{
  if(!TINFOIL){
    double sum_mu0, sum_mu1, sum_mu2;
    sum_mu0 = sum_mu1 = sum_mu2 = 0.0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:sum_mu0, sum_mu1, sum_mu2)
    for(int i = 0; i < nump; i++){
      const double* mui = &mu[i*DIM];
      sum_mu0 += mui[0];
      sum_mu1 += mui[1];
      sum_mu2 += mui[2];
    }

    double dmy_factor = (2.0*M_PI / (2.0*epsilon_bnd + 1.0)) * (cell->iVol);
    energy += dmy_factor*(sum_mu0*sum_mu0 + sum_mu1*sum_mu1 + sum_mu2*sum_mu2);

    dmy_factor *= (-2.0);
#pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < nump; i++){
      int ii = i * DIM;
      const double* mui = &mu[ii];

      efield[ii]   += dmy_factor*sum_mu0;
      efield[ii+1] += dmy_factor*sum_mu1;
      efield[ii+2] += dmy_factor*sum_mu2;
      
      torque[ii]   += dmy_factor*(mui[1]*sum_mu2 - mui[2]*sum_mu1);
      torque[ii+1] += dmy_factor*(mui[2]*sum_mu0 - mui[0]*sum_mu2);
      torque[ii+2] += dmy_factor*(mui[0]*sum_mu1 - mui[1]*sum_mu0);
    }
  }
}
void ewald::compute_r(double &energy, double* force, double* torque, double* efield,
                      double const* r, double const* q, double const* mu,
                      double const* theta) const{

  double dmy_energy = 0.0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy_energy)
  for(int j = 1; j < nump; j++){
    int jj = j*DIM;
    double rij[DIM], cross_mui_muj[DIM], cross_mui_r[DIM], cross_muj_r[DIM];
    double etar, etar2, etar_pi, exp_ewald, erfc_ewald;
    double drij, drij2, drij3, drij5;
    double Br, Cr, Dr;
    double dot_mui_r, dot_muj_r, dot_mui_muj;
    double dmy_force[DIM];

    const double* muj = &mu[jj];
    for(int i = 0; i < j; i++){
      int ii = i*DIM;
      cell->distance_MI(&r[ii], &r[jj], rij);
      drij = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

      if(drij < r2max){
        const double* mui = &mu[ii];
        
        drij = sqrt(drij);
        etar  = eta*drij;
        etar2 = etar*etar;
        etar_pi= 2.0*etar*SQRT_PI_inv;
        exp_ewald = exp(-etar2);
        erfc_ewald = erfc(etar);

        drij = 1.0 / drij;
        drij2= drij*drij;
        drij3= drij2*drij;
        drij5= drij3*drij2;

        Br = (erfc_ewald + etar_pi*exp_ewald) * drij3;
        Cr = (3.0*erfc_ewald + etar_pi*(3.0 + 2.0*etar2)*exp_ewald) * drij5;
        Dr = (15.0*erfc_ewald + etar_pi*(15.0 + 10.0*etar2 + 4.0*etar2*etar2)*exp_ewald) 
	  * drij5 * drij2;

        dot_mui_r = v_inner_prod(mui, rij);
        dot_muj_r = v_inner_prod(muj, rij);
        dot_mui_muj = v_inner_prod(mui, muj);

	v_cross(cross_mui_r, mui, rij);
	v_cross(cross_muj_r, muj, rij);
        v_cross(cross_mui_muj, mui, muj);

        dmy_energy += Br*dot_mui_muj - Cr*dot_mui_r*dot_muj_r;
	for(int d = 0; d < DIM; d++){
	  dmy_force[d] = Cr*(dot_mui_muj*rij[d] + dot_mui_r*muj[d] + dot_muj_r*mui[d])
	    - Dr*(dot_mui_r*dot_muj_r*rij[d]);
	}

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

	//electric fields
#pragma omp atomic
	efield[ii]   += (-Br*muj[0] + Cr*rij[0]*dot_muj_r);
#pragma omp atomic
	efield[ii+1] += (-Br*muj[1] + Cr*rij[1]*dot_muj_r);
#pragma omp atomic
	efield[ii+2] += (-Br*muj[2] + Cr*rij[2]*dot_muj_r);

	efield[jj]   += (-Br*mui[0] + Cr*rij[0]*dot_mui_r);
	efield[jj+1] += (-Br*mui[1] + Cr*rij[1]*dot_mui_r);
	efield[jj+2] += (-Br*mui[2] + Cr*rij[2]*dot_mui_r);
	
	//torques
#pragma omp atomic
	torque[ii]   += (-Br*cross_mui_muj[0] + Cr*cross_mui_r[0]*dot_muj_r);
#pragma omp atomic
	torque[ii+1] += (-Br*cross_mui_muj[1] + Cr*cross_mui_r[1]*dot_muj_r);
#pragma omp atomic
	torque[ii+2] += (-Br*cross_mui_muj[2] + Cr*cross_mui_r[2]*dot_muj_r);

	torque[jj]   += (Br*cross_mui_muj[0] + Cr*cross_muj_r[0]*dot_mui_r);
	torque[jj+1]  += (Br*cross_mui_muj[1] + Cr*cross_muj_r[1]*dot_mui_r);
	torque[jj+2] += (Br*cross_mui_muj[2] + Cr*cross_muj_r[2]*dot_mui_r);
	
      }//rij < r2max
    }//i
  }//j
  energy += dmy_energy;
}

void ewald::precompute_trig_k(double const* r){
  double* iH = (cell->iLambda)[0];
#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < nump; i++){
    const double* ri = &r[i*DIM];

    //transform unit k vectors to cartesian coordinates to compute k.r
    double twopibox_l = PI2*(iH[0]*ri[0] + iH[1]*ri[1]     + iH[2]*ri[2]);
    double twopibox_m = PI2*(iH[DIM]*ri[0] + iH[DIM+1]*ri[1] + iH[DIM+2]*ri[2]);
    double twopibox_n = PI2*(iH[2*DIM]*ri[0] + iH[2*DIM+1]*ri[1] + iH[2*DIM+2]*ri[2]);

    
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
  int sign;

  // cos(p l + q m) = cos(p l) cos(q m) - sin(p l) sin(q m)
  // sin(p l + q m) = sin(p l) cos(q m) + cos(p l) sin(q m)
  sign = (mm > 0 ? 1.0 : -1.0);
  for(int j = 0; j < nump; j++){
    coskr_lm[j] = coskr_l[l][j] * coskr_m[m][j] - sign * sinkr_l[l][j] * sinkr_m[m][j];
    sinkr_lm[j] = sinkr_l[l][j] * coskr_m[m][j] + sign * coskr_l[l][j] * sinkr_m[m][j];
  }

  // cos(p l + q m + r n) = cos(p l + q m) cos(r n) - sin(p l + q m) sin(r n)
  // sin(p l + q m + r n) = sin(p l + q m) cos(r n) + cos(p l + q m) sin(r n)
  sign = (nn > 0 ? 1.0 : -1.0);
  for(int j = 0; j < nump; j++){
    coskr[j] = coskr_lm[j] * coskr_n[n][j] - sign * sinkr_lm[j] * sinkr_n[n][j];
    sinkr[j] = sinkr_lm[j] * coskr_n[n][j] + sign * coskr_lm[j] * sinkr_n[n][j];
  }
}

void ewald::compute_k(double &energy, double* force, double* torque, double* efield, 
                      double const* r, double const* q, double const* mu,
                      double const* theta){
  double dmy_energy, dmy_force, dmy_efield;
  dmy_energy = dmy_force = dmy_efield = 0.0;
  // sum over all k-vectors
  this->precompute_trig_k(r);
  for(int i = 0; i < ewald_domain; i++){
    //compute cos(kr) and sin(kr) terms for all particles
    this->compute_trig_k(ewald_cell[i][0], ewald_cell[i][1], ewald_cell[i][2]);

    double &k0 = ewald_k[i][0];
    double &k1 = ewald_k[i][1];
    double &k2 = ewald_k[i][2];

    double sum_mu_c0, sum_mu_c1, sum_mu_c2;
    double sum_mu_s0, sum_mu_s1, sum_mu_s2;
    sum_mu_c0 = sum_mu_c1 = sum_mu_c2 = 0.0;
    sum_mu_s0 = sum_mu_s1 = sum_mu_s2 = 0.0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:sum_mu_c0, sum_mu_c1, sum_mu_c2, \
							sum_mu_s0, sum_mu_s1, sum_mu_s2)
    for(int j = 0; j < nump; j++){
      int jj = j * DIM;
      const double* pmu = &mu[jj];
      
      sum_mu_c0 += pmu[0] * coskr[j];
      sum_mu_c1 += pmu[1] * coskr[j];
      sum_mu_c2 += pmu[2] * coskr[j];
      
      sum_mu_s0 += pmu[0] * sinkr[j];
      sum_mu_s1 += pmu[1] * sinkr[j];
      sum_mu_s2 += pmu[2] * sinkr[j];
    }

    double sum_k_mu_c = k0 * sum_mu_c0 + k1 * sum_mu_c1 + k2 * sum_mu_c2;
    double sum_k_mu_s = k0 * sum_mu_s0 + k1 * sum_mu_s1 + k2 * sum_mu_s2;
    double kk = k0*k0 + k1*k1 + k2*k2;
    double ewald_damp  = exp(kk*eta_exp)/kk;
    dmy_energy += ewald_damp * (sum_k_mu_c * sum_k_mu_c + sum_k_mu_s * sum_k_mu_s);

    ewald_damp *= (cell->PI8iVol);
#pragma omp parallel for schedule(dynamic, 1) private(dmy_force, dmy_efield)
    for(int j = 0; j < nump; j++){
      int jj = j * DIM;
      const double* pmu = &mu[jj];
      double dot_mu_k = pmu[0]*k0 + pmu[1]*k1 + pmu[2]*k2;

      dmy_force = ewald_damp*dot_mu_k*(sinkr[j]*sum_k_mu_c - coskr[j]*sum_k_mu_s);
      dmy_efield = -ewald_damp*(coskr[j] * sum_k_mu_c + sinkr[j] * sum_k_mu_s);

      force[jj]   += dmy_force * k0;
      force[jj+1] += dmy_force * k1;
      force[jj+2] += dmy_force * k2;

      efield[jj]   += dmy_efield * k0;
      efield[jj+1] += dmy_efield * k1;
      efield[jj+2] += dmy_efield * k2;

      torque[jj]   += dmy_efield * (pmu[1] * k2 - pmu[2] * k1);
      torque[jj+1] += dmy_efield * (pmu[2] * k0 - pmu[0] * k2);
      torque[jj+2] += dmy_efield * (pmu[0] * k1 - pmu[1] * k0);
    }
  }
  energy += (cell->PI4iVol)*dmy_energy;
}

void ewald::compute(double* E_ewald, double* force, double* torque, double* efield,
                    double const* r, double const* q, double const* mu, double const* theta
               ){
  double &e_r = E_ewald[1];
  double &e_k = E_ewald[2];
  double &e_self = E_ewald[3];
  double &e_surface = E_ewald[4];
  e_r = e_k = e_self = e_surface = 0.0;
  clock_t start_t, end_t;
  double cpu_t;

  start_t = clock();
  this->reset(force, torque, efield);
  this->compute_r(e_r, force, torque, efield, r, q, mu, theta);
  this->compute_k(e_k, force, torque, efield, r, q, mu, theta);
  this->compute_self(e_self, force, torque, efield, r, q, mu, theta);
  this->compute_surface(e_surface, force, torque, efield, r, q, mu, theta);
  end_t = clock();
  cpu_t = ((double)end_t - start_t)/CLOCKS_PER_SEC;
  fprintf(stderr, "\tExecution Time = %12.5f \n", cpu_t);
  E_ewald[0] = E_ewald[1] + E_ewald[2] + E_ewald[3] + E_ewald[4];
}

ewald::~ewald(){
  this->free_domain_k();
}
