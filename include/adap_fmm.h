/*
  adap_fmm.h
*/
void adap_fmm_init(const int accuracy, const int nparts);

void adap_fmm_graph(const int nparts, const int s, const double beta, 
		    const double *ploc, const double *pcharge);

void adap_fmm_compute(void);

void adap_fmm_post(double *pot, double *field, double *dxx, double *dyy, double *dzz,
			double *dxy, double *dxz, double *dyz);

void adap_fmm_clean(void);

