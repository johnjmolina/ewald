#ifndef EWALD_H
#define EWALD_H

#include <assert.h>
#include <math.h>
#include <time.h>
#include "parameter_define.h"
#include "alloc.h"
#include "macro.h"
#include "lad3.h"


/*!
  \brief Construct parallelepiped cell reference frame
  \details Basis vector of cell given by cell edges (can be non-uniform)
 */
class parallelepiped{
public:
  parallelepiped(const double a[DIM], const double b[DIM], const double c[DIM]);
  ~parallelepiped();

  /*!
    \brief Compute norm of covariant vector (in cell coordinates)
   */
  double norm_co(const double co_v[DIM]) const;
  /*
    \brief Compute norm of contravariant vector (in cell coordinates)
   */
  double norm_contra(const double contra_v[DIM]) const;
  /*
    \brief Compute minimum image (pbc) distance of vector in lab
    (cartesian) coordinates
   */
  void   distance_MI(const double r1[DIM], const double r2[DIM], 
                     double r12[DIM]) const;
  /*
    \brief Return cell edge lengths
   */
  inline void get_lengths(double &la, double &lb, double &lc) const{
    la = l0;
    lb = l1;
    lc = l2;
  };
  friend class ewald;

private:
  //geometric parameters
  double Vol, iVol;                 //volume, inverse volume
  double PI4iVol, PI8iVol;
  double l0, l1, l2;                //side lengths
  double w0, w1, w2, wmin, wmax;    //perpendicular distances
  double alpha, beta, gamma;        //cosine angles
  
  double Lambda[DIM][DIM];       //transformation matrix
  double tLambda[DIM][DIM];      //transpose
  double iLambda[DIM][DIM];      //inverse
  double tiLambda[DIM][DIM];     //transpose inverse
  double gg[DIM][DIM];           //covariant metric tensor  g_{alpha,beta}
  double GG[DIM][DIM];           //contravariant metric tensor g^{alpha,beta}

  double *e0, *e1, *e2; //basis vectors
  double *E0, *E1, *E2; //dual basis vectors
};

class ewald {
 public:

  /*!
    \brief ewald constructor
   */
  ewald(parallelepiped *_cell,
        const double &ewald_alpha, 
        const double &ewald_epsilon, 
        const double &ewald_delta, 
        const double &ewald_conv,
        const int &num_particles,
        const bool& with_charge, 
        const bool& with_dipole
        );
  /*!
    \brief ewald destroyer
   */
  ~ewald();


  /*!
    \brief Add num_elem particles specified in pid to group with given id
   */
  void add_group(const int &gid, const int &num_elem, const int* pid);

  /*!
    \brief define groups for all particles
   */
  void define_groups(const int* gid);

  /*!
    \brief Reset all particle data
   */
  void reset(double *force, 
	     double *torque, 
	     double *efield,
	     double *efield_grad);

  /*
    \brief Reset boundary (replace current epislon at boundary)
   */
  void reset_boundary(const double &ewald_epsilon);

  /*!
    \brief Compute all particle quantities using multipolar ewald 
   */
  void compute(double* E_ewald,
	       double* force, 
	       double* torque, 
	       double* efield,
	       double* efield_grad,
               double const* r, 
	       double const* q,  
	       double const* mu, 
               const char* save_buffer);

  /*!
    \brief Compute real-space contributions
   */
  void compute_r(double& energy, 
		 double* force, 
		 double* efield,
		 double* efield_grad,
                 double const* r, 
		 double const* q, 
		 double const* mu
		 ) const;

  /*!
    \brief Compute k-space contributions
   */
  void compute_k(double& energy, 
		 double* force, 
		 double* efield,
		 double* efield_grad,
                 double const* r, 
		 double const* q,
		 double const* mu
		 );

  /*!
    \brief Compute surface contributions
   */
  void compute_surface(double &energy, 
		       double* force, 
		       double* efield, 
		       double* efield_grad,
		       double const* r,
		       double const* q, 
		       double const* mu
		       ) const;

  /*!
    \brief Compute self interaction contributions
   */
  void compute_self(double &energy, 
		    double* force, 
		    double* efield, 
		    double* efield_grad,
		    double const* r,
		    double const* q, 
		    double const* mu
		    ) const;

  /*!
    \brief Compute torque on particles given field and field gradient
   */
  void compute_torque(double* torque,
                      double const* efield,
                      double const* efield_grad,
                      double const* mu
		      ) const;

  void compute_upol(double& energy, 
                    double const* polarizability,
                    double const* mu_field) const;

  void compute_mu_induced(double* mu, 
                          double const* polarizability, 
                          double const* efield) const;

 private:
  static const double mu_zero[DIM];

  //particle data
  int nump;
  int* group;

  //ewald parameters
  double rcut, r2max;
  double eta, eta2, eta3, eta_exp;  //screening parameter
  double epsilon_bnd;               //epsilon at boundary
  bool CHARGE, DIPOLE;
  bool TINFOIL;

  // k space parameters
  double k2max;                     //max k norm
  int kmax_l, kmax_m, kmax_n;       //max indices in dual space

  //sekibun cell
  int ewald_domain;                 //number of k vectors
  int** ewald_cell;                 //k vector indices in dual space
  double** ewald_k;                 //k vectors in cartesian coord

  //auxiliary k arrays
  double **coskr_l, **coskr_m, **coskr_n;
  double **sinkr_l, **sinkr_m, **sinkr_n;
  double *coskr_lm, *sinkr_lm;
  double *coskr, *sinkr;

  //cell data
  parallelepiped *cell;

/*!
  \brief Defines the k-space domain that will be used in the ewald calculations
  \details Given the ewald convergence parameters, defines the appropriate cutoff
  for the wave-vectors
  \f[
  k^2 \lessim -4 \eta^2 \delta \ln{\epsilon}
  \f]
  and creates a Sekibun list with all the wave-vectors whose magnitude is less than this cutoff.
  The k vectors are defined as
  \f[
  \vec{k} = 2\pi \left(\bm{h}^{-1}\right)^{\text{t}} \vec{n}
  \f]
  \note For symmetry reasons, we only need to consider half the wave-vectors. We choose the \f$\k_x > 0\f$
  half (for \f$k_x = 0\f$ we take the \f$k_y > 0\f$ part and for \f$k_x=k_y= 0\f$ the \f$k_z > 0\f$ part). 
 */
  void init_domain_k(const double &ewald_delta, const double &ewald_conv);
  void free_domain_k();

  /*!
    \brief Precompute trigonometric sine/cosine factor
   */
  void precompute_trig_k(double const* r);

  /*!
    \brief Compute trigonometric factors for given k vector
   */
  void compute_trig_k(const int&ll, const int&mm, const int&nn);

  /*!
    \brief Compute total charge density in k-space
   */
  void compute_rho_k(double &rho_re, double &rho_im, 
                     double const* q, double const* mu,
                     double const kk[DIM], double const* coskr, double const* sinkr
                     ) const;

  /*!
    \brief auxiliary routines to copy sekibun cell
   */
  void copy_cell(const int &domain, int** old_cell, int** &new_cell);
  void copy_cell(const int &domain, double** old_cell, double** &new_cell);

  /*!
    \brief save data to disk
   */
  void save_results(const double &E_ewald, 
                    double const* force,
                    double const* torque, 
                    double const* efield,
                    double const* efield_grad,
                    double const* r, 
                    double const* q, 
                    double const* mu,
                    char const* save_buffer) const;

  /*!
    \brief save data to disk as cp2k data
   */
  void save_results_cp2k(double const* E_ewald,
                         double const* force,
                         double const* torque,
                         double const* efield,
                         double const* efield_grad,
                         double const* r, 
                         double const* q, 
                         double const* mu,
                         char const* save_buffer) const;
};
#endif
