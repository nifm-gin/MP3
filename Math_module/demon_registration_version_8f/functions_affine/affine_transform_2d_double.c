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

/*
 Affine transformation function (Rotation, Translation, Resize)
 This function transforms a volume with a 3x3 transformation matrix 

 Iout=affine_transform_2d_double(Iin,Minv,mode)

 inputs,
   Iin: The color or greyscale 2D input image
   Minv: The (inverse) 3x3 transformation matrix
   mode: If 0: linear interpolation and outside pixels set to nearest pixel
            1: linear interpolation and outside pixels set to zero
            2: cubic interpolation and outsite pixels set to nearest pixel
            3: cubic interpolation and outside pixels set to zero
 output,
   Iout: The transformed image

 example,
   % Read image
   I=im2double(imread('lenag2.png'))
   % Make a transformation matrix
   M=make_transformation_matrix([2 3],[1.0 1.1],2);
   % Transform the image
   Iout=rigid_transform_2d_double(I,M,1)
   % Show the image
   figure, imshow(Iout);

 Function is written by D.Kroon University of Twente (June 2009)
*/

#ifdef _WIN32
  unsigned __stdcall transformvolume_gray(double **Args) {
#else
  void transformvolume_gray(double **Args) {
#endif
    double *Isize_d, *mean, *A, *Iin, *Iout, *ThreadID, *moded;
    int Isize[3]={0,0,0};
    int mode=0;
    int x,y;
	double *Nthreadsd;
    int Nthreads;

    /* Location of pixel which will be come the current pixel */
    double Tlocalx;
    double Tlocaly;
    
    /* X,Y,Z coordinates of current pixel */
    double xd,yd;

	/* Parts of location calculation */
	double compa0, compa1, compb0, compb1;
	
    /* Variables to store 1D index */
    int indexI;
    
    /* Cubic and outside black booleans */
    bool black, cubic;
    
    /* Multiple threads, one does the odd the other even indexes */
    int offset;
    /* int start; */
    /* int end;; */
    
    Isize_d=Args[0];
    mean=Args[1];
    A=Args[2];
    Iin=Args[3];
    Iout=Args[4];
    ThreadID=Args[5];
    moded=Args[6]; mode=(int) moded[0];
    Nthreadsd=Args[7];  Nthreads=(int)Nthreadsd[0];
   
    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }
    
    Isize[0] = (int)Isize_d[0]; 
    Isize[1] = (int)Isize_d[1]; 
    Isize[2] = (int)Isize_d[2]; 
    
    offset=(int) ThreadID[0];
	
	compb0= A[2] + mean[0];
	compb1= A[5] + mean[1];

    /*  Loop through all image pixel coordinates */
    for (y=offset; y<Isize[1]; y=y+Nthreads)
    {
		yd=(double)y-mean[1];
		compa0 = A[1] *yd + compb0;
		compa1 = A[4] *yd + compb1;

        for (x=0; x<Isize[0]; x++)
        {
            xd=(double)x-mean[0];
            Tlocalx =  A[0] * xd + compa0;
            Tlocaly =  A[3] * xd + compa1;

            /* Set the current pixel value */
            indexI=mindex2(x,y,Isize[0]);
            
            /* interpolate the intensities */
            Iout[indexI]=interpolate_2d_double_gray(Tlocalx, Tlocaly, Isize, Iin,cubic,black); 
        }
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
  unsigned __stdcall transformvolume_color(double **Args) {
#else
  void transformvolume_color(double **Args) {
#endif
    double *Isize_d, *mean, *A, *Iin, *Iout, *ThreadID, *moded;
    int Isize[3]={0,0,0};
    int mode=0;
    int x,y;
	double *Nthreadsd;
    int Nthreads;

    /* Location of pixel which will be come the current pixel */
    double Tlocalx;
    double Tlocaly;
    
    /* X,Y,Z coordinates of current pixel */
    double xd,yd;

	/* Parts of location calculation */
	double compa0, compa1, compb0, compb1;
	
    /* loop throught the colors r,g,b */
    int rgb=0;
            
    /* Variables to store 1D index */
    int indexI;
    
    /* Current voxel/pixel */
    double Ipixel[3]={0,0,0};
    
    /* Cubic and outside black booleans */
    bool black, cubic;
    
    /* Multiple threads, one does the odd the other even indexes */
    int offset;
    /* int start; */
    /* int end;; */
    
    Isize_d=Args[0];
    mean=Args[1];
    A=Args[2];
    Iin=Args[3];
    Iout=Args[4];
    ThreadID=Args[5];
    moded=Args[6]; mode=(int) moded[0];
    Nthreadsd=Args[7];  Nthreads=(int)Nthreadsd[0];
        
    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }
    
    Isize[0] = (int)Isize_d[0]; 
    Isize[1] = (int)Isize_d[1]; 
    Isize[2] = (int)Isize_d[2]; 
    
    offset=(int) ThreadID[0];
	
	compb0= A[2] + mean[0];
	compb1= A[5] + mean[1];

    
    /*  Loop through all image pixel coordinates */
    for (y=offset; y<Isize[1]; y=y+Nthreads)
    {
		yd=(double)y-mean[1];
		compa0 = A[1] *yd + compb0;
		compa1 = A[4] *yd + compb1;

        for (x=0; x<Isize[0]; x++)
        {
            xd=(double)x-mean[0];
            Tlocalx =  A[0] * xd + compa0;
            Tlocaly =  A[3] * xd + compa1;

            /* interpolate the intensities */
            interpolate_2d_double_color(Ipixel,Tlocalx, Tlocaly, Isize, Iin,cubic,black); 
            
            /* Set the current pixel value */
            indexI=mindex2(x,y,Isize[0]);
            for (rgb=0; rgb<3; rgb++)
            {
                Iout[indexI+rgb*Isize[0]*Isize[1]]=Ipixel[rgb];
            }
        }
    }
   
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
                  int nrhs, const mxArray *prhs[] )
{
    /* Ox and Oy are the grid points */
    /* Zo is the input image */
    /* Zi is the transformed image */
    /* nx and ny are the number of grid points (inside the image) */
    double *Iin, *Iout, *M, *moded;
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
    
    /* ID of Threads */
    double **ThreadID;              
    double *ThreadID1;
    
    /* Transformation matrix */
    double A[9]={0,0,0,0,0,0,0,0,0};

    /* Loop variable  */
    int i;
    
    /* Size of input image */
    double Isize_d[3]={0,0,0};
    const mwSize *dims;
    
    double mean[2]={0,0};

  /* Check for proper number of arguments. */
  if(nrhs!=3) {
    mexErrMsgTxt("Three inputs are required.");
  } else if(nlhs!=1) {
    mexErrMsgTxt("One output required");
  }
  
  /* Get the sizes of the image */
  dims = mxGetDimensions(prhs[0]);   
  Isize_d[0] = (double)dims[0]; Isize_d[1] = (double)dims[1]; 
  /* Detect if color image */
  if(mxGetNumberOfDimensions(prhs[0])>2) { Isize_d[2]=(double)3; } else { Isize_d[2]=1; }
  
  /* Create output array */
  if(Isize_d[2]>1) {
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  }
  else  {
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  }
          
  /* Assign pointers to each input. */
  Iin=mxGetPr(prhs[0]);
  M=mxGetPr(prhs[1]);
  moded=mxGetPr(prhs[2]);
 
  A[0] = M[mindex2(0,0,3)]; A[1] = M[mindex2(0,1,3)]; A[2] = M[mindex2(0,2,3)]; 
  A[3] = M[mindex2(1,0,3)]; A[4] = M[mindex2(1,1,3)]; A[5] = M[mindex2(1,2,3)]; 
  A[6] = M[mindex2(2,0,3)]; A[7] = M[mindex2(2,1,3)]; A[8] = M[mindex2(2,2,3)]; 
  
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
	
  
  /* Assign pointer to output. */
  Iout = mxGetPr(plhs[0]);
  
  /* Center of the volume */
  mean[0]=Isize_d[0]/2;  mean[1]=Isize_d[1]/2;  
  
 
  for (i=0; i<Nthreads; i++)
  {
    /*  Make Thread ID  */
    ThreadID1= (double *)malloc( 1* sizeof(double) );
    ThreadID1[0]=i;
    ThreadID[i]=ThreadID1;  
    
	/*  Make Thread Structure  */
    ThreadArgs1 = (double **)malloc( 8* sizeof( double * ) );  
	ThreadArgs1[0]=Isize_d;
	ThreadArgs1[1]=mean;
	ThreadArgs1[2]=A;
	ThreadArgs1[3]=Iin;
	ThreadArgs1[4]=Iout;
	ThreadArgs1[5]=ThreadID[i];
	ThreadArgs1[6]=moded;
	ThreadArgs1[7]=Nthreadsd;
    /* Start a Thread  */
	ThreadArgs[i]=ThreadArgs1;
    if(Isize_d[2]>1) {
		#ifdef _WIN32
			ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume_color, ThreadArgs[i] , 0, NULL );
		#else
			pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &transformvolume_color, ThreadArgs[i]);
		#endif
    }
    else
    {
		#ifdef _WIN32
			ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume_gray, ThreadArgs[i] , 0, NULL );
		#else
			pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &transformvolume_gray, ThreadArgs[i]);
		#endif
    }
  }
   
	#ifdef _WIN32
		for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
		for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
	#else
		for (i=0; i<Nthreads; i++) { pthread_join(ThreadList[i],NULL); }
	#endif


  for (i=0; i<Nthreads; i++) 
  { 
    free(ThreadArgs[i]);
    free(ThreadID[i]);
  }

  free(ThreadArgs);
  free(ThreadID );
  free(ThreadList);

}
        

