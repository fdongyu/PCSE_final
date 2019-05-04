#ifndef functions_h_
#define functions_h_


void init(double x0, double y0, int Npx, int Npy, double *xinit, double *yinit);
double mysecond();
void readtxt(double* xArray, double* yArray, int NR);
void init_random(double x0, double y0, int Npx, int Npy, double *xinit, double *yinit);
int queryUV_index(double xin, double yin, double *xv, double *yv, int NC);

#endif
