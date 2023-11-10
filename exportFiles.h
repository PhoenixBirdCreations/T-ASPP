#include "EOS_cold.h"

void exportEOS_TASPP(char *filename, int n_PT, double *d, double *xm, double *rho_PT, double *a, double *b, double *c, double *f, double extra_cut, double rc, double *Kfitted, double *Gfitted);

void exportQuantities(char *filename, double rc, double fin, struct TASPP eos);
