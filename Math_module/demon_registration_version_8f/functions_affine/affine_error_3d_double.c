#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
/*   undef needed for LCC compiler  */
#undef EXTERN_C
#ifdef _WIN32
  #include <windows.h>
  #include <process.h>
#else 
  #include <pthread.h>
#endif


/* [e,mgrad]=affine_error_3d_double(I1,I2,M,mode);
 * This function transforms a volume with a 4x4 transformation matrix and 
 * also can calculates the error derivatives of the input matrix
 *
 * Function is written by D.Kroon University of Twente (June 2009)
 */

#ifdef _WIN32
  unsigned __stdcall transformvolume_gradient(double **Args) {
#else
  void transformvolume_gradient(double **Args) {
#endif
    double *Isize_d, *mean, *A, *Iin, *Iin2, *ThreadOut, *ThreadID, *moded;
    int Isize[3]={0, 0, 0};
    int x, y, z;
    double *Nthreadsd;
    int Nthreads;
    bool black, cubic;
    int mode=0;
    
    /* Location of pixel which will be come the current pixel */
    double Tlocalx, Tlocaly, Tlocalz;
        
    /* Transformation matrix */
    double B[16]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    double color[8]={0, 0, 0, 0, 0, 0, 0, 0};
    
    /* Total intensity*/
    double I=0;
    
    /* X,Y,Z coordinates of current pixel */
    double xd, yd, zd;
    
    /* Variables to store 1D index */
    int indexI;
    
    /* Delta Finite Difference*/
    double delta=0.00001;
    
    /* Loop variable  */
    int p, i;
    
    /* Multiple threads, one does the odd the other even indexes */
    int offset;
    
    /* Split up matrix multiply to make registration process faster  */
    double acomp0[2]={0,0};
    double acomp1[2]={0,0};
    double acomp2[2]={0,0};
    double bcomp0[2]={0,0};
    double bcomp1[2]={0,0};
    double bcomp2[2]={0,0};
    double ccomp0[2]={0,0};
    double ccomp1[2]={0,0};
    double ccomp2[2]={0,0};
    double dcomp0[2]={0,0}; 
    double dcomp1[2]={0,0}; 
    double dcomp2[2]={0,0};
    
    /* squared difference*/
    double sqd[13]={0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    /* All inputs */
    Isize_d=Args[0];
    mean=Args[1];
    A=Args[2];
    Iin=Args[3];
    Iin2=Args[4];
    ThreadOut=Args[5];
    ThreadID=Args[6];
    Nthreadsd=Args[7];  Nthreads=(int)Nthreadsd[0];
    moded=Args[8]; mode=(int)moded[0];

    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }
    
    
    /* Make gradient matrix */
    for (p=0; p<16; p++)  { B[p]=A[p]+delta; }
    
    Isize[0] = (int)Isize_d[0];
    Isize[1] = (int)Isize_d[1];
    Isize[2] = (int)Isize_d[2];
    
    offset=(int) ThreadID[0];
    
    acomp0[0]=mean[0] + A[3]; acomp1[0]=mean[1] + A[7]; acomp2[0]=mean[2] + A[11];
    acomp0[1]=mean[0] + B[3]; acomp1[1]=mean[1] + B[7]; acomp2[1]=mean[2] + B[11];
    
    /*  Loop through all image pixel coordinates */
    for (z=offset; z<Isize[2]; z=z+Nthreads) 
    {
        zd=z-mean[2];
        bcomp0[0] = A[2] *zd; bcomp1[0] = A[6] *zd; bcomp2[0] = A[10]*zd;
        bcomp0[1] = B[2] *zd; bcomp1[1] = B[6] *zd; bcomp2[1] = B[10]*zd;
        
        for (y=0; y<Isize[1]; y++)
        {
            yd=y-mean[1];
            ccomp0[0] = A[1] *yd; ccomp1[0] = A[5] *yd; ccomp2[0] = A[9] *yd;
            ccomp0[1] = B[1] *yd; ccomp1[1] = B[5] *yd; ccomp2[1] = B[9] *yd; 

            for (x=0; x<Isize[0]; x++) 
            {
                xd=x-mean[0];
                dcomp0[0]=A[0] * xd; dcomp1[0]=A[4] * xd; dcomp2[0]=A[8] * xd;
                dcomp0[1]=B[0] * xd; dcomp1[1]=B[4] * xd; dcomp2[1]=B[8] * xd;
                
                /* Index coordinate current pixel in static image volumeix */
                indexI=mindex3(x, y, z, Isize[0], Isize[1]);
                for(i=0; i<13; i++)
                {
                    Tlocalx = acomp0[0]+bcomp0[0]+ccomp0[0]+dcomp0[0];
                    Tlocaly = acomp1[0]+bcomp1[0]+ccomp1[0]+dcomp1[0];
                    Tlocalz = acomp2[0]+bcomp2[0]+ccomp2[0]+dcomp2[0];
                    switch(i)
                    {
                        case 0: break;
                        case 1: Tlocalx+=dcomp0[1]-dcomp0[0]; break; 
                        case 2: Tlocalx+=ccomp0[1]-ccomp0[0]; break; 
                        case 3: Tlocalx+=bcomp0[1]-bcomp0[0]; break; 
                        case 4: Tlocalx+=acomp0[1]-acomp0[0]; break; 
                        case 5: Tlocaly+=dcomp1[1]-dcomp1[0]; break; 
                        case 6: Tlocaly+=ccomp1[1]-ccomp1[0]; break; 
                        case 7: Tlocaly+=bcomp1[1]-bcomp1[0]; break; 
                        case 8: Tlocaly+=acomp1[1]-acomp1[0]; break; 
                        case 9: Tlocalz+=dcomp2[1]-dcomp2[0]; break; 
                        case 10: Tlocalz+=ccomp2[1]-ccomp2[0]; break; 
                        case 11: Tlocalz+=bcomp2[1]-bcomp2[0]; break;
                        case 12: Tlocalz+=acomp2[1]-acomp2[0]; break;
                    }

                    
                    I=interpolate_3d_double_gray(Tlocalx, Tlocaly, Tlocalz, Isize, Iin, cubic, black);
                    
                    /* Calculate squared difference */
                    sqd[i]+=(I-Iin2[indexI])*(I-Iin2[indexI]);
                }
            }
        }
    }
    
    /* Error and Error gradient to output array */
    for(i=0; i<13; i++)
    {
        ThreadOut[(int)ThreadID[0]*13+i]=sqd[i];
    }
    
    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    #ifdef _WIN32
	_endthreadex( 0 );
    return 0;
	#else
	pthread_exit(NULL);
	#endif
}


#ifdef _WIN32
  unsigned __stdcall transformvolume_error(double **Args) {
#else
  void transformvolume_error(double **Args) {
#endif
    double *Isize_d, *mean, *A, *Iin, *Iin2, *ThreadOut, *ThreadID, *moded;
    int Isize[3]={0, 0, 0};
    int x, y, z;
    double *Nthreadsd;
    int Nthreads;
    bool black, cubic;
    int mode=0;
    
    /* Location of pixel which will be come the current pixel */
    double Tlocalx, Tlocaly, Tlocalz;

      /* Transformation matrix */
    double B[16]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    double color[8]={0, 0, 0, 0, 0, 0, 0, 0};
    
    /* Total intensity*/
    double I=0;
            
    /* X,Y,Z coordinates of current pixel */
    double xd, yd, zd;
    
    /* Variables to store 1D index */
    int indexI;
    
    /* Multiple threads, one does the odd the other even indexes */
    int offset;
    
    /* Split up matrix multiply to make registration process faster  */
    double acomp0, acomp1, acomp2, bcomp0,bcomp1,bcomp2,ccomp0,ccomp1,ccomp2,dcomp0,dcomp1,dcomp2;
    
    /* squared difference */
    double sqd=0;
    
    /* All inputs */
    Isize_d=Args[0];
    mean=Args[1];
    A=Args[2];
    Iin=Args[3];
    Iin2=Args[4];
    ThreadOut=Args[5];
    ThreadID=Args[6];
    Nthreadsd=Args[7];  Nthreads=(int)Nthreadsd[0];
    moded=Args[8]; mode=(int)moded[0];

    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }
       
    Isize[0] = (int)Isize_d[0];
    Isize[1] = (int)Isize_d[1];
    Isize[2] = (int)Isize_d[2];
    
    offset=(int) ThreadID[0];
    
    acomp0=mean[0] + A[3]; 
    acomp1=mean[1] + A[7]; 
    acomp2=mean[2] + A[11];

    /*  Loop through all image pixel coordinates */
    for (z=offset; z<Isize[2]; z=z+Nthreads) 
    {
        zd=z-mean[2];
        bcomp0 = A[2] *zd + acomp0 ; 
        bcomp1 = A[6] *zd + acomp1 ; 
        bcomp2 = A[10]*zd + acomp2 ;
        for (y=0; y<Isize[1]; y++)
        {
            yd=y-mean[1];
            ccomp0 = A[1] *yd + bcomp0; 
            ccomp1 = A[5] *yd + bcomp1; 
            ccomp2 = A[9] *yd + bcomp2;
            for (x=0; x<Isize[0]; x++) 
            {
                xd=x-mean[0];
                dcomp0=A[0] * xd + ccomp0; 
                dcomp1=A[4] * xd + ccomp1; 
                dcomp2=A[8] * xd + ccomp2;
                
                /* Index coordinate current pixel in static image volumeix */
                indexI=mindex3(x, y, z, Isize[0], Isize[1]);

                Tlocalx = dcomp0;
                Tlocaly = dcomp1;
                Tlocalz = dcomp2;
                
                I=interpolate_3d_double_gray(Tlocalx, Tlocaly, Tlocalz, Isize, Iin, cubic, black);
            
                /* Calculate squared difference */
                sqd+=(I-Iin2[indexI])*(I-Iin2[indexI]);
            }
        }
    }
    
    ThreadOut[(int)ThreadID[0]*13]=sqd;
    
    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    #ifdef _WIN32
	_endthreadex( 0 );
    return 0;
	#else
	pthread_exit(NULL);
	#endif
}



/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
    /* Ox and Oy are the grid points */
    /* Zo is the input image */
    /* Zi is the transformed image */
    /* nx and ny are the number of grid points (inside the image) */
    double *Iin, *Iin2, *Ierror, *M, *ThreadOut, *Igradient, *moded;
    mxArray *matlabCallOut[1]={0};
    mxArray *matlabCallIn[1]={0};
    double *Nthreadsd;
    int Nthreads;
    
    /* double pointer array to store all needed function variables) */
    double ***ThreadArgs;
    double **ThreadArgs1;

	/* Handles to the worker threads */
	#ifdef _WIN32
		HANDLE *ThreadList; 
    #else
		pthread_t *ThreadList;
	#endif
	
    
    /* Delta Finite Difference*/
    double delta=0.00001;
    
    /* ID of Threads */
    double **ThreadID;
    double *ThreadID1;
    
    /* Loop variable  */
    int i,j;
    
    /* Transformation matrix */
    double A[16]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    /* Size of input image */
    double Isize_d[3]={0, 0, 0};
    const mwSize *dims;
    
    /*  Size of output */
    int odims[2]={1, 1};
    int odims2[2]={1, 12};
    
    
    double mean[3]={0, 0, 0};
    
    /* Check for proper number of arguments. */
    if(nrhs!=4) {
        mexErrMsgTxt("4 inputs are required.");
    } 

    if(nlhs==2) /* Also gradient needed */
    {
        /* Create output array */
        plhs[1] = mxCreateNumericArray(2, odims2, mxDOUBLE_CLASS, mxREAL);
        /* Connect output array */
        Igradient = (double *)mxGetData(plhs[1]);
    }
     
    /* nsubs=mxGetNumberOfDimensions(prhs[0]);  */
    
    /* Get the sizes of the image */
    dims = mxGetDimensions(prhs[0]);
    Isize_d[0] = (double)dims[0]; Isize_d[1] = (double)dims[1]; Isize_d[2] = (double)dims[2];
    
    
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    
    /* Assign pointers to each input. */
    Iin=(double *)mxGetData(prhs[0]);
    Iin2=(double *)mxGetData(prhs[1]);
    M=(double *)mxGetData(prhs[2]);
    moded=(double *)mxGetData(prhs[3]);
            
    A[0] = M[mindex2(0, 0, 4)];  A[1] = M[mindex2(0, 1, 4)];  A[2] = M[mindex2(0, 2, 4)];  A[3] = M[mindex2(0, 3, 4)];
    A[4] = M[mindex2(1, 0, 4)];  A[5] = M[mindex2(1, 1, 4)];  A[6] = M[mindex2(1, 2, 4)];  A[7] = M[mindex2(1, 3, 4)];
    A[8] = M[mindex2(2, 0, 4)];  A[9] = M[mindex2(2, 1, 4)];  A[10] = M[mindex2(2, 2, 4)]; A[11] = M[mindex2(2, 3, 4)];
    A[12] = M[mindex2(3, 0, 4)]; A[13] = M[mindex2(3, 1, 4)]; A[14] = M[mindex2(3, 2, 4)]; A[15] = M[mindex2(3, 3, 4)];
    
    mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
    Nthreadsd=mxGetPr(matlabCallOut[0]);
    Nthreads=(int)Nthreadsd[0];
    /* Reserve room for handles of threads in ThreadList  */
	#ifdef _WIN32
		ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
    #else
		ThreadList = (pthread_t*)malloc(Nthreads* sizeof( pthread_t ));
	#endif
	
	
    ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
    ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );
    ThreadOut = (double *)malloc( 13*Nthreads* sizeof(double) );
    
    /* Assign pointer to output. */
    Ierror = (double *)mxGetData(plhs[0]);
    
    /* Center of the volume */
    mean[0]=Isize_d[0]/2;  mean[1]=Isize_d[1]/2;  mean[2]=Isize_d[2]/2;
    
    for (i=0; i<Nthreads; i++) {
        /*  Make Thread ID  */
        ThreadID1= (double *)malloc( 1* sizeof(double) );
        ThreadID1[0]=(double)i;
        ThreadID[i]=ThreadID1;
        
        /*  Make Thread Structure  */
        ThreadArgs1 = (double **)malloc(9* sizeof( double * ) );
        ThreadArgs1[0]=Isize_d;
        ThreadArgs1[1]=mean;
        ThreadArgs1[2]=A;
        ThreadArgs1[3]=Iin;
        ThreadArgs1[4]=Iin2;
        ThreadArgs1[5]=ThreadOut;
        ThreadArgs1[6]=ThreadID[i];
        ThreadArgs1[7]=Nthreadsd;
        ThreadArgs1[8]=moded;
        
        /* Start a Thread  */
        ThreadArgs[i]=ThreadArgs1;
        if(nlhs==2)
        {
			#ifdef _WIN32
					ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume_gradient, ThreadArgs[i] , 0, NULL );
			#else
					pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &transformvolume_gradient, ThreadArgs[i]);
			#endif
	
        }
        else
        {
			#ifdef _WIN32
					ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume_error, ThreadArgs[i] , 0, NULL );
			#else
					pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &transformvolume_error, ThreadArgs[i]);
			#endif

        }
    }
    
	
	#ifdef _WIN32
		for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
		for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
	#else
		for (i=0; i<Nthreads; i++) { pthread_join(ThreadList[i],NULL); }
	#endif
    
	Ierror[0]=0;
    for (i=0; i<Nthreads; i++)
    {
        Ierror[0]+=ThreadOut[i*13];
    }
    Ierror[0]/= Isize_d[0]*Isize_d[1]*Isize_d[2];
    if(nlhs==2) /* Also gradient needed */
    {
        for (j=0; j<12; j++)
        {
			Igradient[j]=0;
            for (i=0; i<Nthreads; i++) 
            {
                Igradient[j]+=ThreadOut[i*13+j+1]-ThreadOut[i*13];
            }
            Igradient[j]/= Isize_d[0]*Isize_d[1]*Isize_d[2]*delta;
        }
    }
    
    
    for (i=0; i<Nthreads; i++) {
        free(ThreadArgs[i]);
        free(ThreadID[i]);
    }
    
    free(ThreadArgs);
    free(ThreadOut);
    free(ThreadID );
    free(ThreadList);
}


