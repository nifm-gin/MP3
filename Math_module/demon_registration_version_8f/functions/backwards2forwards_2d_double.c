#include "mex.h"
#include "math.h"
/*
 This function will turn a backwards transformation field into 
 a forwards transformation field.

 [Fx,Fy]=backwards2forwards_2d_double(Bx,By,H); 

 inputs,
   Bx,By : The backward transformation fields
   H : The Splatting kernel 
 outputs,
   Fx,Fy : The forward transformation fields

 Function is written by D.Kroon University of Twente (Februari 2009)
*/


/* Convert 2D/3D matrix index to 1D index */
int mindex2(int x, int y, int sizx) { return y*sizx+x; }

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *Fx, *Fy, *Bx, *By, *H, *Num;
    
    /* Size of input transformation fields */
    mwSize  Bsizex, Bsizey;
    const mwSize *Bdims;

    /* Size of the input kernel */
    mwSize  Hsizex, Hsizey;
    const mwSize *Hdims;
    int hx_center,  hy_center;
     
    /* Variables to store 1D index */
    int index;

    /* Variable to store vector */
    double valx, valy;
    
    /* Loop variables */
    int i, t, iHx, iHy;
        
    /* Linear interpolation variables */
    int xBas0, xBas1,yBas0,yBas1;
    double perc[4], perct;
    double xCom, yCom;
    double Tlocalx, Tlocaly;
    
    /* X,Y coordinates of current pixel */
    int x,y, tx, ty;
    
    /* Check for proper number of arguments. */
    if(nrhs!=3) {
        mexErrMsgTxt("Three inputs are required.");
    } else if(nlhs!=2) {
        mexErrMsgTxt("Two outputs are required");
    }
  
  /* Assign pointers to each input. */
  Bx=mxGetPr(prhs[0]);
  By=mxGetPr(prhs[1]);
  H=mxGetPr(prhs[2]);
 
  /* Get the sizes of the kernel  */
  Hdims = mxGetDimensions(prhs[2]);   
  Hsizex = Hdims[0]; Hsizey = Hdims[1];
  
  /* Get the sizes of the input transformation fields */
  Bdims = mxGetDimensions(prhs[0]);   
  Bsizex = Bdims[0]; Bsizey = Bdims[1];

  /* Create array to count number of kernels added */
  Num = (double *)malloc(Bsizex*Bsizey*sizeof(double));  
  for (i=0; i<(Bsizex*Bsizey); i++){ Num[i] = 0;}    

  /* Create output arrays */
  plhs[0] = mxCreateNumericArray(2, Bdims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(2, Bdims, mxDOUBLE_CLASS, mxREAL);

  /* Assign pointer to output. */
  Fx = mxGetPr(plhs[0]);
  Fy = mxGetPr(plhs[1]);

  
  /* Gaussian kernel center */
  hx_center=(int)-floor((double)Hsizex/2); 
  hy_center=(int)-floor((double)Hsizey/2); 
  
  /* Loop through all image pixel coordinates */
  for (y=0; y<Bsizey; y++)
  {
    for (x=0; x<Bsizex; x++)
    {
        valx=-Bx[mindex2(x,y,Bsizex)];
        valy=-By[mindex2(x,y,Bsizex)];
        Tlocalx =x-valx;
        Tlocaly =y-valy;
   
        /* Determine the coordinates of the pixel(s) which will be come the current pixel */
        /* (using linear interpolation)   */
        xBas0=(int) floor(Tlocalx); 
        yBas0=(int) floor(Tlocaly);
        xBas1=xBas0+1;      
        yBas1=yBas0+1;

        /* Linear interpolation constants (percentages) */
        xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);
        perc[0]=(1-xCom) * (1-yCom);
        perc[1]=(1-xCom) * yCom;
        perc[2]=xCom * (1-yCom);
        perc[3]=xCom * yCom;
        
        /* Loop through the whole kernel */
        for (iHy=0; iHy<Hsizey; iHy++)
        {
            for (iHx=0; iHx<Hsizex; iHx++)
            {
                /* Process all 4 neighbors                 */
                for(t=0; t<4; t++)
                {
                   if(t==0){
                        tx=xBas0+iHx+hx_center; ty=yBas0+iHy+hy_center; 
                        perct=perc[0]*H[mindex2(iHx,iHy,Hsizex)];
                   }
                   else if(t==1){
                        tx=xBas0+iHx+hx_center; ty=yBas1+iHy+hy_center; 
                        perct=perc[1]*H[mindex2(iHx,iHy,Hsizex)];
                   }
                   else if(t==2){
                        tx=xBas1+iHx+hx_center; ty=yBas0+iHy+hy_center; 
                        perct=perc[2]*H[mindex2(iHx,iHy,Hsizex)];
                   }
                   else{
                        tx=xBas1+iHx+hx_center; ty=yBas1+iHy+hy_center; 
                        perct=perc[3]*H[mindex2(iHx,iHy,Hsizex)];
                   }
                   if((tx>=0)&&(ty>=0)&&(tx<Bsizex)&&(ty<Bsizey))
                   {
                        index=mindex2(tx,ty,Bsizex);
                        Fx[index]+=valx*perct; 
                        Fy[index]+=valy*perct;
                        Num[index]=Num[index]+perct;
                   }
               }
                
            }
        }
    }
  }
  for (y=0; y<Bsizey; y++)
  {
      for (x=0; x<Bsizex; x++)
      {
        index=mindex2(x,y,Bsizex);
        Fx[index]/=(Num[index]+0.00000001); 
        Fy[index]/=(Num[index]+0.00000001);
      }
  }
  free(Num);
}
        

