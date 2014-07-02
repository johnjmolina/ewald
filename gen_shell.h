#ifndef GEN_SHELL_H
#define GEN_SHELL_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <assert.h>
#include "macro.h"

using namespace std;

bool rcompare(vector <int> a, vector <int> b);

//Generate list of cencentric spherical cells 
void nshell_list(const int &ncut, const double &la, const double &lb, const double &lc,
		 int &lmax, int &mmax, int &nmax, long int &ntot,
		 vector<double> &nmag, vector< vector< vector <int> > > &nlist,
                 FILE* fout);
#endif
