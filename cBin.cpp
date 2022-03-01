//-------------------------------------------------------------------------------
// Bins vector X into N equal compartments.
//
//-------------------------------------------------------------------------------


#include <math.h>
#include "mex.h"

#include <stdio.h>
#include <iostream>
/* input */
#define IN	prhs[0]
#define Nr	prhs[1]
#define Nc	prhs[2]

/* output */
#define OUT  plhs[0]

// sum data
double sum(double * x, int first, int last){
    double s = 0;
    for(int i = first; i <= last; i ++)
        s += x[i];
    return s;
}

// bin 1D image
void bin(double * x, int len, double * out, int N){
    
    double *splits;
    splits = new double [N+1];
    for(int j=0;j<N+1;j++)
            splits[j] = 1+(double)len/N*j;
    
    for(int k=0;k<N;k++){
        if(k == 0)
            out[k] = sum(x,(int)splits[k]-1,floor(splits[k+1])-2) + x[(int)floor(splits[k+1])-1]*(splits[k+1] - floor(splits[k+1]));
        else if(k == N-1)
            out[k] = sum(x,floor(splits[k])-1,splits[k+1]-2) - x[(int)floor(splits[k])-1]*(splits[k] - floor(splits[k]));
        else 
            out[k] = sum(x,floor(splits[k])-1,floor(splits[k+1])-2) + x[(int)floor(splits[k+1])-1]*(splits[k+1] - floor(splits[k+1]))- x[(int)floor(splits[k])-1]*(splits[k] - floor(splits[k]));
    }

  delete [] splits;
}

//-------------------------------------------------------------------------------
// bin 2D image
void bin2(double * in, int sr, int sc, double * out, int nr, int nc){
    
    double * interm;
    interm = new double [nr*sc];
    double * raw_r, * raw_c, * binned_r, * binned_c;
    raw_r = new double [sr];
    raw_c = new double [sc];
    binned_r = new double [nr];
    binned_c = new double [nc];

    
    
    for(int i=0;i<sc;i++){
        for(int j=0;j<sr;j++)
            raw_r[j] = in[i*sr+j];
        bin(raw_r,sr,binned_r,nr);
        for(int j=0;j<nr;j++)
            interm[i*nr+j] = binned_r[j];
    }
    
    
    for(int i=0;i<nr;i++){
        for(int j=0;j<sc;j++)
            raw_c[j] = interm[j*nr+i];
        bin(raw_c,sc,binned_c,nc);
        for(int j=0;j<nc;j++)
            out[j*nr+i] = binned_c[j];
    }
    
    
    delete [] raw_r;
    delete [] raw_c;
    delete [] binned_r;
    delete [] binned_c;
    delete [] interm;
}


/* mex function */
//-------------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[], 
          int nrhs, const mxArray*prhs[] )
{ 
    double * in; 
    double * out; 
    int sr, sc; 
    
    /* Check format */    
    if (nrhs != 3) { 
    mexErrMsgTxt("only three inputs"); 
	} 
    else if (nlhs > 1) {
    mexErrMsgTxt("only one output"); 
    } 

    /* calculate image size */ 
    sr = mxGetM(IN); 
    sc = mxGetN(IN);


    /* get pointer */ 
    in = mxGetPr(IN);
    int nr       = (double)(mxGetScalar(Nr));
    int nc       = (double)(mxGetScalar(Nc));
    
    /* create matrix */ 
    OUT = mxCreateDoubleMatrix(nr, nc, mxREAL); 
    out = mxGetPr(OUT);
    
    /* bin data */
	bin2(in, sr, sc, out, nr, nc);

 //   return;
}

