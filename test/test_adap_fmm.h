/*
  test_adap_fmm.h
*/

void test_parser(int argc, char **argv, double *beta, int *nparts, int *s, 
		 int *accuracy, int *distribution);

void test_init(const int distribution, const int nparts, const int accuracy, 
	       const int s, const double beta, 
	       double **ploc, double **pcharge, double **pot, double **field,
	       double **dxx, double **dyy, double **dzz, double **dxy, 
	       double **dxz, double **dyz);

void test_verify(const double beta, const int nparts, const double *ploc, 
		 const double *pcharge, const double *pot, const double *field, 
		 const int accuracy, const double *dxx, const double *dyy, 
		 const double *dzz, const double *dxy, const double *dxz,
		 const double *dyz);

void test_clean(double *ploc, double *pcharge, double *pot, double *field, double *dxx,
		double *dyy, double *dzz, double *dxy, double *dxz, double *dyz);

void test_data(const int nparts, const int distribution, double *ploc, double *pcharge);

void lapdirect(const int nparts, const double *ploc, const double *pcharge, 
	       const int pid, double *pot, double *field, double *dxx, 
	       double *dyy, double *dzz, double *dxy, double *dxz, double *dyz);

void yukdirect(const double beta, const int nparts, const double *ploc, 
	       const double *pcharge, const int i, double *pot, double *field,
	       double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz);

#define NDIRECT 400

