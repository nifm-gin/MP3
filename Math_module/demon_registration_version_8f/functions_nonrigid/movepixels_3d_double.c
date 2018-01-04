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

/*  This function movepixels, will translate the pixels of an image
 *  according to x, y and z translation images (bilinear interpolated). 
 * 
 *  Iout = movepixels_3d_double(I,Tx,Ty,Tz);
 *
 *  Function is written by D.Kroon University of Twente (July 2009)
 */

#ifdef _WIN32
  unsigned __stdcall transformvolume(double **Args) {
#else
  void transformvolume(double **Args) {
#endif
    /* I is the input image, Iout the transformed image  */
    /* Tx and Ty images of the translation of every pixel. */
    double *Iin, *Iout, *Tx, *Ty, *Tz;
    double *Nthreadsd;
    int Nthreads;
	/*  if one outside pixels are set to zero. */
	double  *moded;
	int mode=0;
   
    /* Cubic and outside black booleans */
    bool black, cubic;
    
    /* 3D index storage*/
    int indexI;
    
    /* Size of input image */
    double *Isize_d;
    int Isize[3]={0,0,0};
        
    /* Location of translated pixel */
    double Tlocalx;
    double Tlocaly;
    double Tlocalz;
    
    /* offset */
    int ThreadOffset=0;
    
    /* The thread ID number/name */
    double *ThreadID;
    
    /* X,Y coordinates of current pixel */
    int x,y,z;
    
    Iin=Args[0];
    Iout=Args[1];
    Tx=Args[2];
    Ty=Args[3];
    Tz=Args[4];
    Isize_d=Args[5];
    ThreadID=Args[6];
	moded=Args[7]; mode=(int) moded[0];
	Nthreadsd=Args[8];  Nthreads=(int)Nthreadsd[0];

    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }
    
    Isize[0] = (int)Isize_d[0]; 
    Isize[1] = (int)Isize_d[1]; 
    Isize[2] = (int)Isize_d[2]; 
       
    ThreadOffset=(int) ThreadID[0];
	
    /*  Loop through all image pixel coordinates */
    for (z=ThreadOffset; z<Isize[2]; z=z+Nthreads)
	{
        for (y=0; y<Isize[1]; y++)
        {
            for (x=0; x<Isize[0]; x++)
            {
                Tlocalx=((double)x)+Tx[mindex3(x,y,z,Isize[0],Isize[1])];
                Tlocaly=((double)y)+Ty[mindex3(x,y,z,Isize[0],Isize[1])];
                Tlocalz=((double)z)+Tz[mindex3(x,y,z,Isize[0],Isize[1])];
                
                /* Set the current pixel value */
                indexI=mindex3(x,y,z,Isize[0],Isize[1]);
                Iout[indexI]=interpolate_3d_double_gray(Tlocalx, Tlocaly, Tlocalz, Isize, Iin,cubic,black); 
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
    /* I is the input image, Iout the transformed image  */
    /* Tx and Ty images of the translation of every pixel. */
    double *Iin, *Iout, *Tx, *Ty, *Tz;
	double *moded;
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

    /* Size of input image */
    mwSize Isizex, Isizey,Isizez;
    double Isize_d[3]={0,0,0};
    const mwSize *dims;

	/* Loop variable  */
	int i;
	
    /* Check for proper number of arguments. */
    if(nrhs!=5) {
      mexErrMsgTxt("Five inputs are required.");
    } else if(nlhs!=1) {
      mexErrMsgTxt("One output required");
    }
 
    /* Get the sizes of the input image */
    dims = mxGetDimensions(prhs[0]);   
    Isizex = (mwSize)dims[0]; 
    Isizey = (mwSize)dims[1];
    Isizez = (mwSize)dims[2];

    Isize_d[0]=Isizex;  Isize_d[1]=Isizey; Isize_d[2]=Isizez;
    
    /* Create image matrix for the return arguments with the size of input image  */  
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); 
    
    /* Assign pointers to each input. */
    Iin=(double *)mxGetData(prhs[0]);
    Tx=(double *)mxGetData(prhs[1]);
    Ty=(double *)mxGetData(prhs[2]);
    Tz=(double *)mxGetData(prhs[3]);
	moded=(double *)mxGetData(prhs[4]);
	

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
    Iout = (double *)mxGetData(plhs[0]);
   
  for (i=0; i<Nthreads; i++)
  {
    /*  Make Thread ID  */
    ThreadID1= (double *)malloc( 1* sizeof(double) );
    ThreadID1[0]=i;
    ThreadID[i]=ThreadID1;  
    
	/*  Make Thread Structure  */
    ThreadArgs1 = (double **)malloc( 9* sizeof( double * ) );  
    ThreadArgs1[0]=Iin;
    ThreadArgs1[1]=Iout;
    ThreadArgs1[2]=Tx;
    ThreadArgs1[3]=Ty;
    ThreadArgs1[4]=Tz;
    ThreadArgs1[5]=Isize_d;
    ThreadArgs1[6]=ThreadID[i];
	ThreadArgs1[7]=moded;
	ThreadArgs1[8]=Nthreadsd;
	
    /* Start a Thread  */
	ThreadArgs[i]=ThreadArgs1;
	#ifdef _WIN32
		ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs[i] , 0, NULL );
	#else
		pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &transformvolume, ThreadArgs[i]);
	#endif
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
        

