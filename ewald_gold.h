#ifndef DIPOLE_GOLD_H
#define DIPOLE_GOLD_H
#include <iostream>
#include <vector>
#include <assert.h>
#include <time.h>
#include "macro.h"
#include "alloc.h"
#include "parameter_define.h"
#include "gen_shell.h"
#include "ewald.h"


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
                             );

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
                       );
#endif
