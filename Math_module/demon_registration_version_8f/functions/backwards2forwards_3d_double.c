#include "mex.h"
#include "math.h"
/*
 This function will turn a backwards transformation field into 
 a forwards transformation field.

 [Fx,Fy,Fz]=backwards2forwards_2d_double(Bx,By,Bz,H); 

 inputs,
   Bx,By,Bz : The backward transformation fields
   H : The Splatting kernel 
 outputs,
   Fx,Fy,Fz : The forward transformation fields

 Function is written by D.Kroon University of Twente (Februari 2009)
*/


/* Convert 2D/3D matrix index to 1D index */
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *Fx, *Fy, *Fz, *Bx, *By, *Bz, *H, *Num;
    
    /* Size of input transformation fields */
    mwSize  Bsizex, Bsizey, Bsizez;
    const mwSize *Bdims;

    /* Size of the input kernel */
    mwSize  Hsizex, Hsizey, Hsizez;
    const mwSize *Hdims;
    int hx_center,  hy_center, hz_center;
     
    /* Variables to store 1D index */
    int index;

    /* Variable to store vector */
    double valx, valy, valz;
    
    /* Loop variables */
    int i, t, iHx, iHy, iHz;
        
    /* Linear interpolation variables */
    int xBas0, xBas1,yBas0,yBas1,zBas0,zBas1;
    double perc[8], perct;
    double xCom, yCom, zCom;
    double Tlocalx, Tlocaly, Tlocalz;
    
    /* X,Y coordinates of current pixel */
    int x,y, z, tx, ty, tz;
    
    /* Check for proper number of arguments. */
    if(nrhs!=4) {
        mexErrMsgTxt("Four inputs are required.");
    } else if(nlhs!=3) {
        mexErrMsgTxt("Three outputs are required");
    }
  
  /* Get the sizes of the kernel  */
  Hdims = mxGetDimensions(prhs[3]);   
  Hsizex = Hdims[0]; Hsizey = Hdims[1]; Hsizez = Hdims[2];
  
  /* Get the sizes of the input transformation fields */
  Bdims = mxGetDimensions(prhs[0]);   
  Bsizex = Bdims[0]; Bsizey = Bdims[1]; Bsizez = Bdims[2];
    
  /* Create output arrays */
  plhs[0] = mxCreateNumericArray(3, Bdims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(3, Bdims, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericArray(3, Bdims, mxDOUBLE_CLASS, mxREAL);

  /* Create array to count number of kernels added */
  Num = (double *)malloc(Bsizex*Bsizey*Bsizez*sizeof(double));  
  for (i=0; i<(Bsizex*Bsizey*Bsizez); i++){ Num[i] = 0;}    
  /* Assign pointers to each input. */
  Bx=mxGetPr(prhs[0]);
  By=mxGetPr(prhs[1]);
  Bz=mxGetPr(prhs[2]);
  H=mxGetPr(prhs[3]);
 
  /* Assign pointer to output. */
  Fx = mxGetPr(plhs[0]);
  Fy = mxGetPr(plhs[1]);
  Fz = mxGetPr(plhs[2]);
    
  /* Gaussian kernel center */
  hx_center=(int)-floor((double)Hsizex/2); 
  hy_center=(int)-floor((double)Hsizey/2); 
  hz_center=(int)-floor((double)Hsizez/2); 
  
  /* Loop through all image pixel coordinates */
  for (z=0; z<Bsizez; z++)
  {
      for (y=0; y<Bsizey; y++)
      {
        for (x=0; x<Bsizex; x++)
        {
            valx=-Bx[mindex3(x,y,z,Bsizex,Bsizey)];
            valy=-By[mindex3(x,y,z,Bsizex,Bsizey)];
            valz=-Bz[mindex3(x,y,z,Bsizex,Bsizey)];
            
            Tlocalx =x-valx;
            Tlocaly =y-valy;
            Tlocalz =z-valz;
            
            /* Determine the coordinates of the pixel(s) which will be come the current pixel */
            /* (using linear interpolation)   */
            xBas0=(int) floor(Tlocalx); 
            yBas0=(int) floor(Tlocaly);
            zBas0=(int) floor(Tlocalz);
            xBas1=xBas0+1;      
            yBas1=yBas0+1;
            zBas1=zBas0+1;

            /* Linear interpolation constants (percentages) */
            xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);  zCom=Tlocalz-floor(Tlocalz);
            perc[0]=(1-xCom) * (1-yCom) * (1-zCom);
            perc[1]=(1-xCom) * (1-yCom) * zCom;
            perc[2]=(1-xCom) * yCom * (1-zCom);
            perc[3]=(1-xCom) * yCom * zCom;
            perc[4]=xCom * (1-yCom) * (1-zCom);
            perc[5]=xCom * (1-yCom) * zCom;
            perc[6]=xCom * yCom * (1-zCom);
            perc[7]=xCom * yCom * zCom;

            /* Loop through the whole kernel */
            for (iHz=0; iHz<Hsizez; iHz++)
            {
                for (iHy=0; iHy<Hsizey; iHy++)
                {
                    for (iHx=0; iHx<Hsizex; iHx++)
                    {
                       /* Process all 4 neighbors                 */
                        for(t=0; t<8; t++)
                        {
                           if(t==0){
                                tx=xBas0+iHx+hx_center; ty=yBas0+iHy+hy_center; tz=zBas0+iHz+hz_center; 
                                perct=perc[0]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }
                           else if(t==1){
                                tx=xBas0+iHx+hx_center; ty=yBas0+iHy+hy_center; tz=zBas1+iHz+hz_center; 
                                perct=perc[1]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }
                           else if(t==2){
                                tx=xBas0+iHx+hx_center; ty=yBas1+iHy+hy_center; tz=zBas0+iHz+hz_center; 
                                perct=perc[2]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }
                           else if(t==3){
                                tx=xBas0+iHx+hx_center; ty=yBas1+iHy+hy_center; tz=zBas1+iHz+hz_center; 
                                perct=perc[3]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }
                           else if(t==4){
                                tx=xBas1+iHx+hx_center; ty=yBas0+iHy+hy_center; tz=zBas0+iHz+hz_center; 
                                perct=perc[4]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }
                           else if(t==5){
                                tx=xBas1+iHx+hx_center; ty=yBas0+iHy+hy_center; tz=zBas1+iHz+hz_center; 
                                perct=perc[5]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }
                           else if(t==6){
                                tx=xBas1+iHx+hx_center; ty=yBas1+iHy+hy_center; tz=zBas0+iHz+hz_center; 
                                perct=perc[6]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }
                           else{
                                tx=xBas1+iHx+hx_center; ty=yBas1+iHy+hy_center; tz=zBas1+iHz+hz_center; 
                                perct=perc[7]*H[mindex3(iHx,iHy,iHz,Hsizex,Hsizey)];
                           }

                           if((tx>=0)&&(ty>=0)&&(tz>=0)&&(tx<Bsizex)&&(ty<Bsizey)&&(tz<Bsizez))
                           {
                                index=mindex3(tx,ty,tz,Bsizex,Bsizey);
                                Fx[index]+=valx*perct; 
                                Fy[index]+=valy*perct;
                                Fz[index]+=valz*perct;
                                Num[index]=Num[index]+perct;
                           }
                        }
                   }
                }
            }
        }
      }
  }
  for (z=0; z<Bsizez; z++)
  {
      for (y=0; y<Bsizey; y++)
      {
          for (x=0; x<Bsizex; x++)
          {
            index=mindex3(x,y,z,Bsizex,Bsizey);
            Fx[index]/=(Num[index]+0.00000001); 
            Fy[index]/=(Num[index]+0.00000001);
            Fz[index]/=(Num[index]+0.00000001);
          }
      }
  }
  free(Num);
}
        

