#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define pi 3.1415926
//#include "functions.h"

void init(double x0, double y0, int Npx, int Npy, double *xinit, double *yinit){
    
    int i, j;
    double *x=malloc(Npx*sizeof(double)), *y=malloc(Npy*sizeof(double));
    double xlen=100, ylen=100; // unit: m
    double x_interval, y_interval;
    
    x[0]=x0-xlen/2.;
    y[0]=y0-ylen/2.;
    x_interval=xlen/Npx;
    y_interval=ylen/Npy;
    
    for (i=1;i<Npx;i++){
        x[i]=x[0]+x_interval*i;
    }
    for (i=1;i<Npy;i++){
        y[i]=y[0]+y_interval*i;
    }
    //printf("%6.13f\n", x0);
    for (i=0;i<Npx;i++){
        for (j=0;j<Npy;j++){
            xinit[i*Npy+j] = x[i];
            yinit[i*Npy+j] = y[j];
            //printf("%6.13f %6.13f\n", xinit[i*10+j], yinit[i*10+j]);
        }
    }

    printf("%6.13f %6.13f\n", xinit[55], yinit[55]);

}


void init_random(double x0, double y0, int Npx, int Npy, double *xinit, double *yinit){

    double R = Npx*Npy; //radius
    int N = Npx*Npy;
    double r, t;
    
    int i;
    for (i=0;i<N;i++){
        t = 2*pi*random() / (float)RAND_MAX;
        r = R * sqrt(random() / (float)RAND_MAX);
        xinit[i] = x0+r*cos(t);
        yinit[i] = y0+r*sin(t);
    }
    
}




double CalculateDistance(x1, y1, x2, y2){
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

int find_minimum(double *a, int n) {
    int c, min, index;
    min = a[0];
    index = 0;
    
    for (c = 1; c < n; c++) {
        if (a[c] < min) {
            index = c;
            min = a[c];
        }
    }
    //printf("index is %d \n", (int)index);
    return index;
}


int queryUV_index(double xin, double yin, double *xv, double *yv, int NC){
    
    double *distance=malloc(NC*sizeof(double));
    int i;
    for (i=0;i<NC;i++){
        distance[i]=CalculateDistance(xv[i],yv[i],xin,yin);
    }
    return find_minimum(distance, NC);
}


double mysecond()
{
    struct timeval tp;
    struct timezone tzp;
    int i;
    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


void readtxt(double* xArray, double* yArray, int NR){
    
    FILE *myFile;
    myFile = fopen("../txt/cells.dat", "r");
    
    //read file into array
    //int NR = 1000; // release locations
    int Nall = NR * 8;
    double  *temArray=malloc(Nall*sizeof(double));
    //double  *xArray=malloc(NR*sizeof(double));
    //double  *yArray=malloc(NR*sizeof(double));
    int i,j;
    
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);
    }
    
    for (i = 0; i < Nall; i++){
        fscanf(myFile, "%lf", &temArray[i]);
    }
    
    for (i = 0; i< NR; i++){
        j = i*8;
        xArray[i] = temArray[j];
        yArray[i] = temArray[j+1];
        //printf("Number is: %6.5f, %6.5f\n", xArray[i], yArray[i]);
        }
    
    fclose(myFile);
}
