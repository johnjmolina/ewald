#include "gen_shell.h"

bool rcompare(vector <int> a, vector <int> b){
  assert(a.size() == 3 && b.size() == 3);
  double r2a = (double) (a.at(0)*a.at(0) + a.at(1)*a.at(1) + a.at(2)*a.at(2));
  double r2b = (double) (b.at(0)*b.at(0) + b.at(1)*b.at(1) + b.at(2)*b.at(2));
  return (r2a < r2b);
}

//Generate list of cencentric spherical cells 
void nshell_list(const int &ncut, const double &la, const double &lb, const double &lc,
		 int &lmax, int &mmax, int &nmax, long int &ntot,
		 vector<double> &nmag, vector< vector< vector <int> > > &nlist,
                 FILE* fout){
  
  double da, db, dc, dl, dr;
  double diag = sqrt(la*la + lb*lb + lc*lc); // length of diagonal
  double dlmax = (double)ncut * diag;
  int dmy_lmax = (int)(dlmax/la) + 2;
  int dmy_mmax = (int)(dlmax/lb) + 2;
  int dmy_nmax = (int)(dlmax/lc) + 2;
  lmax = 0;
  mmax = 0;
  nmax = 0;

  int nshells = ncut + 1;
  nmag.resize(nshells);
  nlist.resize(nshells);
  fprintf(fout, "#Generate spherical image cell list for direct space calculations \n");
  fprintf(fout, "Number of shells: %d\n", nshells);
  fprintf(fout, "lcut (dr): %f (%f)\n", dlmax, diag);
  for(int i = 0; i < nshells; i++){
    dl = (double)i * diag;
    nmag.at(i) = dl; 
  }

  for(int ll = -dmy_lmax; ll <= dmy_lmax; ll++){
    da = (double)ll * la;

    for(int mm = -dmy_mmax; mm <= dmy_mmax; mm++){
      db = (double)mm * lb;

      for(int nn = -dmy_nmax; nn <= dmy_nmax; nn++){
	dc = (double)nn * lc;
	dr = sqrt(da*da + db*db + dc*dc);

	if(dr <= dlmax){
	  lmax = MAX(lmax, abs(ll));
	  mmax = MAX(mmax, abs(mm));
	  nmax = MAX(nmax, abs(nn));

	  int dmy_n[] = {ll, mm, nn};
	  vector<int> dmy_vector(dmy_n, dmy_n + sizeof(dmy_n)/sizeof(int));

	  int i = 0;
	  while(dr > nmag.at(i))
	    i++;
	  if(nlist.at(i).size() == 0){
	    vector< vector<int> > dmy_list;
	    dmy_list.push_back(dmy_vector);
	    nlist.at(i) = dmy_list;
	  }else{
	    nlist.at(i).push_back(dmy_vector);
	  }

	} //dr < dlmax

      }
    }
  }
  fprintf(fout, " lmax : %d\n", lmax);
  fprintf(fout, " mmax : %d\n", mmax);
  fprintf(fout, " nmax : %d\n", nmax);
  fprintf(fout, "\n");
  
  //sort intermediate shell arrays
  for(int i = 0; i < nshells; i++){
    vector< vector<int> > dmy_list(nlist.at(i));
    sort(dmy_list.begin(), dmy_list.end(), rcompare);
    nlist.at(i) = dmy_list;
  }

  //print results
  ntot = 0;
  for(int i = 0; i < nshells; i++){
    ntot += (long int) nlist.at(i).size();
  }
  fprintf(fout, "Total number of cells: %ld\n\n", ntot);
}
