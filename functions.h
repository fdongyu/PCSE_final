#ifndef functions_h_
#define functions_h_


void init(float x0, float y0, int Npx, int Npy, float *xinit, float *yinit);
double mysecond();
void readtxt(float* xArray, float* yArray, int NR);
void init_random(float x0, float y0, int Npx, int Npy, float *xinit, float *yinit);
int queryUV_index(float xin, float yin, float *xv, float *yv, int NC);

#endif
