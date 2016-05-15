/*
  RPY.h
*/

void test_parser(int argc, char ** argv, double *beta, int *nparts, int *s, 
		 int *accuracy, int *distribution);

void test_init(const int distribution, const int nparts, const int accuracy, 
	       const int s, const double beta, double **ploc, double **pcharge, 
	       double **pot, double **field, double **dxx, double **dyy, double **dzz,
	       double **dxy, double **dxz, double **dyz);

void test_clean(double *ploc, double *pcharge, double *charge, double *pot, double *RPY, double *field,
		  double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz );

void test_data(int nparts, int distribution, double *ploc, double *pcharge);


