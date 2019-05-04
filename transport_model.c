#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>

#include "functions.h"

#define NDIMS 2
#define NDIMS1 1

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#define NC 21484
#define NT 97
#define NK 20


double *offset(double *array, int i, int j, int width) {
    return &array[width * i + j];
}

void copyArray(double* array1, double* array2, int n){
    
    int i;
    //printf("Number of Particles:  %d \n",n);
    for (i=0;i<n;i++){
        array1[i] = array2[i];
    }
}



void main(){
    
    /* Timing information */
    double i0, i1, i2, i3, i4, i5;
    double t0, t1, t2, t3;
    
    i0 = mysecond();
    
    /* Read u,v,t from netcdf file */
    int ncid_in, varid;
    int retval;
    double *xv=malloc(NC*sizeof(double)), *yv=malloc(NC*sizeof(double));
    
    static size_t start1d[] = {0}, start3d[] = {0, 0, 0}; /* start at first value */
    static size_t count1d[] = {NC}, count3d[] = {NT, NK, NC};
    static size_t count1d_nt[] = {NT};
    static double outuc[(int)NT][(int)NK][(int)NC], outvc[(int)NT][(int)NK][(int)NC];
    double *time=malloc(NT*sizeof(double));
    char FILE_NAME[] = "DATA/GalvCoarse_0000.nc";
    
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid_in)))
    ERR(retval);
    
    /* Read the time.  */
    if ((retval = nc_inq_varid(ncid_in, "time", &varid)))
    ERR(retval);
    if ((retval = nc_get_vara_double(ncid_in, varid, start1d, count1d_nt, &time[0])))
    ERR(retval);
    
    /* Read the x.  */
    if ((retval = nc_inq_varid(ncid_in, "xv", &varid)))
    ERR(retval);
    if ((retval = nc_get_vara_double(ncid_in, varid, start1d, count1d, &xv[0])))
    ERR(retval);
    
    /* Read the y. */
    if ((retval = nc_inq_varid(ncid_in, "yv", &varid)))
    ERR(retval);
    if ((retval = nc_get_vara_double(ncid_in, varid, start1d, count1d, &yv[0])))
    ERR(retval);
    
    /* Read the uc. */
    if ((retval = nc_inq_varid(ncid_in, "uc", &varid)))
    ERR(retval);
    if ((retval = nc_get_vara_double(ncid_in, varid, start3d, count3d, &outuc[0][0][0])))
    ERR(retval);
    
    /* Read the vc. */
    if ((retval = nc_inq_varid(ncid_in, "vc", &varid)))
    ERR(retval);
    if ((retval = nc_get_vara_double(ncid_in, varid, start3d, count3d, &outvc[0][0][0])))
    ERR(retval);
    
    printf("uc is %6.13f\n", (double)outuc[88][0][21000]);
    printf("vc is %6.13f\n", (double)outvc[88][0][21000]);
    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid_in)))
    ERR(retval);
    
    printf("*** SUCCESS reading example file %s!\n", FILE_NAME);
    
    
    i1 = mysecond();
    
    /* transport model core */
    
    /* Step 1: define particle initial locations*/
    /* 10*10 matrix */
    int Npx=10, Npy=10;
    int NP=Npx*Npy; // number of particles at one location
    double *xtem=malloc(NP*sizeof(double)), *ytem=malloc(NP*sizeof(double));
    double x0, y0;
    int NR=5; // number of release locations
    int ir; // release location iterater
    int Ninit = NR*NP;
    double *xinit=malloc(Ninit*sizeof(double)), *yinit=malloc(Ninit*sizeof(double)); // row:NR, col: NP
    
    double  *xArray=malloc(NR*sizeof(double));
    double  *yArray=malloc(NR*sizeof(double));
    readtxt(xArray, yArray, NR);
    
    /*
    int ii;
    for (ii = 0; ii< NR; ii++){
        printf("Number is: %6.5f, %6.5f\n", xArray[ii], yArray[ii]);
    }
    */
    
    for (ir = 0; ir< NR; ir++){
        x0 = xArray[ir];
        y0 = yArray[ir];
        init_random(x0, y0, Npx, Npy, xtem, ytem);
        copyArray(offset(xinit, ir, 0, NP), xtem, NP);
        copyArray(offset(yinit, ir, 0, NP), ytem, NP);
    }

    free(xtem);
    free(ytem);
    free(xArray);
    free(yArray);
    
    //init(x0, y0, Npx, Npy, xinit, yinit);
    
    i2 = mysecond();
    
    /* Step 2: particle transport calculation */
    int n, k, j, i;
    int period=48; // transport period, units: hours
    double dt=720; // unit: seconds
    double utem, vtem, dx, dy;
    int Ntp=3600/dt*period; // Ntp: particle time length
    int ind0, ind1, ind2, ind3;
    int tind0, tind1, tind2, tind3;
    double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;
    double *time_p=malloc(Ntp*sizeof(double)); //time for particle transport
    double x[Ninit][Ntp], y[Ninit][Ntp];
    //double *x=malloc(Ninit*Ntp*sizeof(double)), *y=malloc(Ninit*Ntp*sizeof(double)); //particle  trajectories
    
    i3 = mysecond();
    
    for (i=0;i<Ninit;i++){
        x[i][0]=xinit[i];
        y[i][0]=yinit[i];
    }
    
    time_p[0] = time[1] - 3600;
    
    for (i=0;i<Ninit;i++){
        printf("Particle %d \n",(int)i);
        for (n=1;n<Ntp;n++){
            /* 1st step */
            ind0=queryUV_index(x[i][n-1],y[i][n-1], xv, yv, NC);
            tind0=(int)n/(3600/dt);
            //printf("%6.13f  %6.13f; index is %d, tindex is %d\n", x[n-1], y[n-1], (int)ind, (int)tind);
            utem=outuc[tind0][0][ind0]; vtem=outvc[tind0][0][ind0];
            
            
            /*
            k1x=dt*utem; k1y=dt*vtem;
            
            ind1=queryUV_index(x[i][n-1]+k1x/2.,y[i][n-1]+k1y/2., xv, yv, NC);
            tind1=(int)(n+0.5)/(3600/dt);
            utem=outuc[tind1][0][ind1]; vtem=outvc[tind1][0][ind1];
            k2x=dt*utem; k2y=dt*vtem;
            
            ind2=queryUV_index(x[i][n-1]+k2x/2.,y[i][n-1]+k2y/2., xv, yv, NC);
            tind2=(int)(n+0.5)/(3600/dt);
            utem=outuc[tind2][0][ind2]; vtem=outvc[tind2][0][ind2];
            k3x=dt*utem; k3y=dt*vtem;
            
            ind3=queryUV_index(x[i][n-1]+k3x,y[i][n-1]+k3y, xv, yv, NC);
            tind3=(int)(n+1.)/(3600/dt);
            utem=outuc[tind3][0][ind3]; vtem=outvc[tind3][0][ind3];
            k4x=dt*utem; k4y=dt*vtem;
            
            dx = (k1x+2*k2x+2*k3x+k4x)/6.; dy = (k1y+2*k2y+2*k3y+k4y)/6.;
            */
            
            dx=dt*utem; dy=dt*vtem;
            
            
            x[i][n] = x[i][n-1] + dx;
            y[i][n] = y[i][n-1] + dy;
            time_p[n] = time_p[n-1] + dt;
            //printf("x is %6.13f; y is %6.13f\n", (double)x[n], (double)y[n]);
        }
    }
    
    i4 = mysecond();
    
    FILE *stream;
    stream=fopen("DATA/init.dat","w");
    for (n=0;n<Ninit;n++){
        fprintf(stream,"%6.13f  %6.13f\n",xinit[n], yinit[n]);}
    fclose(stream);
    
    
    /* write to netcdf file*/
    //int retval;
    int ncid_out, x_dimid, y_dimid, np_dimid, nt_dimid, varid_x, varid_y, varid_t;
    int dimids[NDIMS], dimids1[NDIMS1];
    size_t chunks[NDIMS];
    int shuffle, deflate, deflate_level;
    /* Set chunking, shuffle, and deflate. */
    shuffle = NC_SHUFFLE;
    deflate = 1;
    deflate_level = 1;
    
    /* Create the file. The NC_NETCDF4 flag tells netCDF to
     * create a netCDF-4/HDF5 file.*/
    if ((retval = nc_create("DATA/xy_Euler.nc", NC_NETCDF4, &ncid_out)))
    ERR(retval);
    
    /* Define the dimensions in the root group. Dimensions are visible
     * in all subgroups. */
    if ((retval = nc_def_dim(ncid_out, "nt", Ntp, &nt_dimid)))
    ERR(retval);
    if ((retval = nc_def_dim(ncid_out, "np", Ninit, &np_dimid)))
    ERR(retval);
    
    
    dimids[0] = np_dimid;
    chunks[0] = Ninit/4;
    dimids[1] = nt_dimid;
    chunks[1] = Ntp/4;
    
    dimids1[0] = nt_dimid;
    
    
    /* Define the variable. */
    if ((retval = nc_def_var(ncid_out, "X", NC_DOUBLE, NDIMS,
                             dimids, &varid_x)))
    ERR(retval);
    
    
    /*
    if ((retval = nc_def_var_chunking(ncid_out, varid_x, 0, &chunks[0])))
    ERR(retval);
    if ((retval = nc_def_var_deflate(ncid_out, varid_x, shuffle, deflate,
                                     deflate_level)))
    ERR(retval);
    */
    
    if ((retval = nc_def_var(ncid_out, "Y", NC_DOUBLE, NDIMS,
                             dimids, &varid_y)))
    ERR(retval);
    
    if ((retval = nc_def_var(ncid_out, "time", NC_DOUBLE, NDIMS1,
                             dimids1, &varid_t)))
    ERR(retval);
    
    /* No need to explicitly end define mode for netCDF-4 files. Write
     * the pretend data to the file. */
    if ((retval = nc_put_var_double(ncid_out, varid_x, &x[0][0])))
    ERR(retval);
    if ((retval = nc_put_var_double(ncid_out, varid_y, &y[0][0])))
    ERR(retval);
    if ((retval = nc_put_var_double(ncid_out, varid_t, &time_p[0])))
    ERR(retval);
    
    /* Close the file. */
    if ((retval = nc_close(ncid_out)))
    ERR(retval);
    
    printf("*** SUCCESS writing example file xy.nc!\n");
    
    i5 = mysecond();
    
    /* timing */
    t0 = i1-i0;   // Read input file
    t1 = i2-i1;   // Initialization
    t2 = i4-i3;   // particle tracking
    t3 = i5-i4;   // Save output file

    
    printf("Action          ::     time/s     Time resolution = 1.0E-04\n");
    printf("-------\n");
    printf("CPU: Read Input File    ::     %13.3f (sec)\n", t0);
    printf("CPU: Initialization     ::     %13.3f (sec)\n", t1);
    printf("CPU: Particle Tracking  ::     %13.3f (sec)\n", t2);
    printf("CPU: Save Output file   ::     %13.3f (sec)\n", t3);
    
}


