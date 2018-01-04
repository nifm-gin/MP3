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
 *  according to x and y translation images (bilinear interpolated). 
 * 
 *  Iout = movepixels_2d_double(I,Tx,Ty);
 *
 *  Function is written by D.Kroon University of Twente (July 2009)
 */

#ifdef _WIN32
  unsigned __stdcall transformvolume_color(double **Args) {
#else
  void transformvolume_color(double **Args) {
#endif

    /* I is the input image, Iout the transformed image  */
    /* Tx and Ty images of the translation of every pixel. */
    double *Iin, *Iout, *Tx, *Ty;
    double *Nthreadsd;
    int Nthreads;
	/*  if one outside pixels are set to zero. */
	double  *moded;
	int mode=0;
	
    /* 2D index storage */
    int indexI;
    
    /* Size of input image */
    int Isize[3]={0,0,0};
    double *Isize_d;
    
    /* Location of translated pixel */
    double Tlocalx;
    double Tlocaly;
    
    /* Cubic and outside black booleans */
    bool black, cubic;

    /* loop throught the colors r,g,b */
    int rgb=0;
    
    /* Current voxel/pixel */
    double Ipixel[3]={0,0,0};
    
    /* offset */
    int offset=0;
    
    /* The thread ID number/name */
    double *ThreadID;
    
    /* X,Y coordinates of current pixel */
    int x,y;
    
    Iin=Args[0];
    Iout=Args[1];
    Tx=Args[2];
    Ty=Args[3];
    Isize_d=Args[4];
    ThreadID=Args[5];
	moded=Args[6]; mode=(int) moded[0];
    Nthreadsd=Args[7];  Nthreads=(int)Nthreadsd[0];
	
    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }

    Isize[0]=(int)Isize_d[0];
    Isize[1]=(int)Isize_d[1];
    Isize[2]=(int)Isize_d[2];
        
    offset=(int) ThreadID[0];
	
    /*  Loop through all image pixel coordinates */
    for (y=offset; y<Isize[1]; y=y+Nthreads)
    {
        for (x=0; x<Isize[0]; x++)
        {
            Tlocalx=((double)x)+Tx[mindex2(x,y,Isize[0])];
            Tlocaly=((double)y)+Ty[mindex2(x,y,Isize[0])];
            
            /* interpolate the intensities */
            interpolate_2d_double_color(Ipixel,Tlocalx, Tlocaly, Isize, Iin,cubic,mode);
         
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


#ifdef _WIN32
  unsigned __stdcall transformvolume_gray(double **Args) {
#else
  void transformvolume_gray(double **Args) {
#endif
    /* I is the input image, Iout the transformed image  */
    /* Tx and Ty images of the translation of every pixel. */
    double *Iin, *Iout, *Tx, *Ty;
    double *Nthreadsd;
    int Nthreads;
	/*  if one outside pixels are set to zero. */
	double  *moded;
	int mode=0;

   
    /* Cubic and outside black booleans */
    bool black, cubic;

    /* 2D index storage */
    int indexI;
    
    /* Size of input image */
    int Isize[3]={0,0,0};
    double *Isize_d;
    
    /* Location of translated pixel */
    double Tlocalx;
    double Tlocaly;
    
    /* loop throught the colors r,g,b */
    int rgb=0;
    
    /* Current voxel/pixel */
    double Ipixel[3]={0,0,0};
    
    /* offset */
    int offset=0;
    
    /* The thread ID number/name */
    double *ThreadID;
    
    /* X,Y coordinates of current pixel */
    int x,y;
    
    Iin=Args[0];
    Iout=Args[1];
    Tx=Args[2];
    Ty=Args[3];
    Isize_d=Args[4];
    ThreadID=Args[5];
	moded=Args[6]; mode=(int) moded[0];
    Nthreadsd=Args[7];  Nthreads=(int)Nthreadsd[0];
	  
    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }

    Isize[0]=(int)Isize_d[0];
    Isize[1]=(int)Isize_d[1];
    Isize[2]=(int)Isize_d[2];
        
    offset=(int) ThreadID[0];
	
    /*  Loop through all image pixel coordinates */
    for (y=offset; y<Isize[1]; y=y+Nthreads)
    {
        for (x=0; x<Isize[0]; x++)
        {
            Tlocalx=((double)x)+Tx[mindex2(x,y,Isize[0])];
            Tlocaly=((double)y)+Ty[mindex2(x,y,Isize[0])];
            
            /* Set the current pixel value */
            indexI=mindex2(x,y,Isize[0]);
           
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


/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* I is the input image, Iout the transformed image  */
    /* Tx and Ty images of the translation of every pixel. */
    double *Iin, *Iout, *Tx, *Ty;
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
    const mwSize *dims;
    double Isize_d[3]={0,0,0};
    
	/* Loop variable  */
	int i;
	
    /* Check for proper number of arguments. */
    if(nrhs!=4) {
      mexErrMsgTxt("Four inputs are required.");
    } else if(nlhs!=1) {
      mexErrMsgTxt("One output required");
    }
 
    /* Get the sizes of the image */
    dims = mxGetDimensions(prhs[0]);   
    Isize_d[0] = (double)dims[0]; Isize_d[1] = (double)dims[1]; 
    /* Detect if color image */
    if(mxGetNumberOfDimensions(prhs[0])>2) { Isize_d[2]=(double)3; } else { Isize_d[2]=1; }
  
    /* Create image matrix for the return arguments with the size of input image  */  
    if(Isize_d[2]>1) {
          plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    }
    else  {
        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    }

    /* Assign pointers to each input. */
    Iin=mxGetPr(prhs[0]);
    Tx=mxGetPr(prhs[1]);
    Ty=mxGetPr(prhs[2]);
	moded=mxGetPr(prhs[3]);
	  

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
   
  for (i=0; i<Nthreads; i++)
  {
    /*  Make Thread ID  */
    ThreadID1= (double *)malloc( 1* sizeof(double) );
    ThreadID1[0]=i;
    ThreadID[i]=ThreadID1;  
    
	/*  Make Thread Structure  */
    ThreadArgs1 = (double **)malloc( 8* sizeof( double * ) );  
	ThreadArgs1[0]=Iin;
    ThreadArgs1[1]=Iout;
    ThreadArgs1[2]=Tx;
    ThreadArgs1[3]=Ty;
    ThreadArgs1[4]=Isize_d;
    ThreadArgs1[5]=ThreadID[i];
	ThreadArgs1[6]=moded;
	ThreadArgs1[7]=Nthreadsd;
	
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
        

