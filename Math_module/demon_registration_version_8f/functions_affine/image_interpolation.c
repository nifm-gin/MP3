#include "math.h"
#include "image_interpolation.h"

/* Image and Volume interpolation 
 *
 * Function is written by D.Kroon University of Twente (June 2009)
 */
 

/* Get an pixel from an image, if outside image, black or nearest pixel */
double getintensity_mindex2(int x, int y, int sizx, int sizy, double *I) {
    return I[y*sizx+x];
}

/* Get an pixel from an image, if outside image, black or nearest pixel */
double getcolor_mindex2(int x, int y, int sizx, int sizy, double *I, int rgb) {
    return I[rgb*sizy*sizx+y*sizx+x];
}

/* Get an pixel from an image, if outside image, black or nearest pixel */
double getcolor_mindex3(int x, int y, int z, int sizx, int sizy, int sizz, double *I) {
    return I[z*sizx*sizy+y*sizx+x];
}

/* Get an pixel from an image, if outside image, black or nearest pixel */
float getcolor_mindex3_float(int x, int y, int z, int sizx, int sizy, int sizz, float *I) {
    return I[z*sizx*sizy+y*sizx+x];
}

double interpolate_2d_linear_gray_black(double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1;
    double perc[4]={0, 0, 0, 0};
    double xCom, yCom, xComi, yComi;
    double Ipixel=0;
    
    
    /*  Rounded location  */
    double fTlocalx, fTlocaly;
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    xBas1=xBas0+1; yBas1=yBas0+1;
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx; yCom=Tlocaly-fTlocaly;
    xComi=(1-xCom); yComi=(1-yCom);
    perc[0]=xComi * yComi; perc[1]=xComi * yCom; perc[2]=xCom * yComi; perc[3]=xCom * yCom;
    
    if((xBas0>=0)&&(xBas0<Isize[0])) {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            Ipixel+=getintensity_mindex2(xBas0, yBas0, Isize[0], Isize[1], Iin)*perc[0];
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            Ipixel+=getintensity_mindex2(xBas0, yBas1, Isize[0], Isize[1], Iin)*perc[1];
        }
    }
    if((xBas1>=0)&&(xBas1<Isize[0]))  {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            Ipixel+=getintensity_mindex2(xBas1, yBas0, Isize[0], Isize[1], Iin)*perc[2];
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            Ipixel+=getintensity_mindex2(xBas1, yBas1, Isize[0], Isize[1], Iin)*perc[3];
        }
    }
    return Ipixel;
}

double interpolate_2d_cubic_gray_black(double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /* Floor of coordinate */
    double fTlocalx, fTlocaly;
    /* Zero neighbor */
    int xBas0, yBas0;
    /* The location in between the pixels 0..1 */
    double tx, ty;
    /* Neighbor loccations */
    int xn[4], yn[4];
    
    /* The vectors */
    double vector_tx[4], vector_ty[4];
    double vector_qx[4], vector_qy[4];
    /* Interpolated Intensity; */
    double Ipixel=0, Ipixelx=0;
    /* Loop variable */
    int i;
    
    /* Determine of the zero neighbor */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly;
    
    /* Determine the t vectors */
    vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
    vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
    vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
    vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
    vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
    vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
    vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
    vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
    vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for(i=0; i<4; i++) {
        Ipixelx=0;
        if((yn[i]>=0)&&(yn[i]<Isize[1])) {
            if((xn[0]>=0)&&(xn[0]<Isize[0])) {
                Ipixelx+=vector_qx[0]*getintensity_mindex2(xn[0], yn[i], Isize[0], Isize[1], Iin);
            }
            if((xn[1]>=0)&&(xn[1]<Isize[0])) {
                Ipixelx+=vector_qx[1]*getintensity_mindex2(xn[1], yn[i], Isize[0], Isize[1], Iin);
            }
            if((xn[2]>=0)&&(xn[2]<Isize[0])) {
                Ipixelx+=vector_qx[2]*getintensity_mindex2(xn[2], yn[i], Isize[0], Isize[1], Iin);
            }
            if((xn[3]>=0)&&(xn[3]<Isize[0])) {
                Ipixelx+=vector_qx[3]*getintensity_mindex2(xn[3], yn[i], Isize[0], Isize[1], Iin);
            }
        }
        Ipixel+= vector_qy[i]*Ipixelx;
    }
    return Ipixel;
}

double interpolate_2d_cubic_gray(double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /* Floor of coordinate */
    double fTlocalx, fTlocaly;
    /* Zero neighbor */
    int xBas0, yBas0;
    /* The location in between the pixels 0..1 */
    double tx, ty;
    /* Neighbor loccations */
    int xn[4], yn[4];
    
    /* The vectors */
    double vector_tx[4] , vector_ty[4];
    double vector_qx[4], vector_qy[4];
    /* Interpolated Intensity; */
    double Ipixel=0, Ipixelx=0;
    /* Temporary value boundary */
    int b;
    /* Loop variable */
    int i;
    
    /* Determine of the zero neighbor */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly;
    
    /* Determine the t vectors */
    vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
    vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
    vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
    vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
    vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
    vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
    vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
    vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
    vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    
    /* Clamp to image boundary if outside image */
    if(xn[0]<0) { xn[0]=0;if(xn[1]<0) { xn[1]=0;if(xn[2]<0) { xn[2]=0; if(xn[3]<0) { xn[3]=0; }}}}
    if(yn[0]<0) { yn[0]=0;if(yn[1]<0) { yn[1]=0;if(yn[2]<0) { yn[2]=0; if(yn[3]<0) { yn[3]=0; }}}}
    b=Isize[0]-1;
    if(xn[3]>b) { xn[3]=b;if(xn[2]>b) { xn[2]=b;if(xn[1]>b) { xn[1]=b; if(xn[0]>b) { xn[0]=b; }}}}
    b=Isize[1]-1;
    if(yn[3]>b) { yn[3]=b;if(yn[2]>b) { yn[2]=b;if(yn[1]>b) { yn[1]=b; if(yn[0]>b) { yn[0]=b; }}}}
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for(i=0; i<4; i++) {
        Ipixelx =vector_qx[0]*getintensity_mindex2(xn[0], yn[i], Isize[0], Isize[1], Iin);
        Ipixelx+=vector_qx[1]*getintensity_mindex2(xn[1], yn[i], Isize[0], Isize[1], Iin);
        Ipixelx+=vector_qx[2]*getintensity_mindex2(xn[2], yn[i], Isize[0], Isize[1], Iin);
        Ipixelx+=vector_qx[3]*getintensity_mindex2(xn[3], yn[i], Isize[0], Isize[1], Iin);
        Ipixel+= vector_qy[i]*Ipixelx;
    }
    return Ipixel;
}

double interpolate_2d_linear_gray(double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1;
    double perc[4]={0, 0, 0, 0};
    double xCom, yCom, xComi, yComi;
    double color[4]={0, 0, 0, 0};
    
    /*  Rounded location  */
    double fTlocalx, fTlocaly;
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    xBas1=xBas0+1; yBas1=yBas0+1;
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx; yCom=Tlocaly-fTlocaly;
    xComi=(1-xCom); yComi=(1-yCom);
    perc[0]=xComi * yComi;
    perc[1]=xComi * yCom;
    perc[2]=xCom * yComi;
    perc[3]=xCom * yCom;
    
    if(xBas0<0) { xBas0=0; if(xBas1<0) { xBas1=0; }}
    if(yBas0<0) { yBas0=0; if(yBas1<0) { yBas1=0; }}
    if(xBas1>(Isize[0]-1)) { xBas1=Isize[0]-1; if(xBas0>(Isize[0]-1)) { xBas0=Isize[0]-1; }}
    if(yBas1>(Isize[1]-1)) { yBas1=Isize[1]-1; if(yBas0>(Isize[1]-1)) { yBas0=Isize[1]-1; }}
    
    color[0]=getintensity_mindex2(xBas0, yBas0, Isize[0], Isize[1], Iin);
    color[1]=getintensity_mindex2(xBas0, yBas1, Isize[0], Isize[1], Iin);
    color[2]=getintensity_mindex2(xBas1, yBas0, Isize[0], Isize[1], Iin);
    color[3]=getintensity_mindex2(xBas1, yBas1, Isize[0], Isize[1], Iin);
    return color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
}

void interpolate_2d_cubic_color_black(double *Ipixel, double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /* Floor of coordinate */
    double fTlocalx, fTlocaly;
    /* Zero neighbor */
    int xBas0, yBas0;
    /* The location in between the pixels 0..1 */
    double tx, ty;
    /* Neighbor loccations */
    int xn[4], yn[4];
    
    /* The vectors */
    double vector_tx[4], vector_ty[4];
    double vector_qx[4], vector_qy[4];
    double Ipixelx[3];
    /* Loop variable */
    int i;
    
    /* Determine of the zero neighbor */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly;
    
    /* Determine the t vectors */
    vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
    vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
    vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
    vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
    vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
    vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
    vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
    vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
    vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    Ipixel[0]=0; Ipixel[1]=0; Ipixel[2]=0;
    for(i=0; i<4; i++) {
        Ipixelx[0]=0; Ipixelx[1]=0; Ipixelx[2]=0;
        if((yn[i]>=0)&&(yn[i]<Isize[1])) {
            if((xn[0]>=0)&&(xn[0]<Isize[0])) {
                Ipixelx[0]+=vector_qx[0]*getcolor_mindex2(xn[0], yn[i], Isize[0], Isize[1], Iin, 0);
                Ipixelx[1]+=vector_qx[0]*getcolor_mindex2(xn[0], yn[i], Isize[0], Isize[1], Iin, 1);
                Ipixelx[2]+=vector_qx[0]*getcolor_mindex2(xn[0], yn[i], Isize[0], Isize[1], Iin, 2);
            }
            if((xn[1]>=0)&&(xn[1]<Isize[0])) {
                Ipixelx[0]+=vector_qx[1]*getcolor_mindex2(xn[1], yn[i], Isize[0], Isize[1], Iin, 0);
                Ipixelx[1]+=vector_qx[1]*getcolor_mindex2(xn[1], yn[i], Isize[0], Isize[1], Iin, 1);
                Ipixelx[2]+=vector_qx[1]*getcolor_mindex2(xn[1], yn[i], Isize[0], Isize[1], Iin, 2);
            }
            if((xn[2]>=0)&&(xn[2]<Isize[0])) {
                Ipixelx[0]+=vector_qx[2]*getcolor_mindex2(xn[2], yn[i], Isize[0], Isize[1], Iin, 0);
                Ipixelx[1]+=vector_qx[2]*getcolor_mindex2(xn[2], yn[i], Isize[0], Isize[1], Iin, 1);
                Ipixelx[2]+=vector_qx[2]*getcolor_mindex2(xn[2], yn[i], Isize[0], Isize[1], Iin, 2);
            }
            if((xn[3]>=0)&&(xn[3]<Isize[0])) {
                Ipixelx[0]+=vector_qx[3]*getcolor_mindex2(xn[3], yn[i], Isize[0], Isize[1], Iin, 0);
                Ipixelx[1]+=vector_qx[3]*getcolor_mindex2(xn[3], yn[i], Isize[0], Isize[1], Iin, 1);
                Ipixelx[2]+=vector_qx[3]*getcolor_mindex2(xn[3], yn[i], Isize[0], Isize[1], Iin, 2);
            }
        }
        Ipixel[0]+= vector_qy[i]*Ipixelx[0];
        Ipixel[1]+= vector_qy[i]*Ipixelx[1];
        Ipixel[2]+= vector_qy[i]*Ipixelx[2];
    }
}

void interpolate_2d_cubic_color(double *Ipixel, double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /* Floor of coordinate */
    double fTlocalx, fTlocaly;
    /* Zero neighbor */
    int xBas0, yBas0;
    /* The location in between the pixels 0..1 */
    double tx, ty;
    /* Neighbor loccations */
    int xn[4], yn[4];
    
    /* The vectors */
    double vector_tx[4], vector_ty[4];
    double vector_qx[4], vector_qy[4];
    double Ipixelx;
    /* Temporary value boundary */
    int b;
    /* Loop variable */
    int i, rgb;
    
    /* Determine of the zero neighbor */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly;
    
    /* Determine the t vectors */
    vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
    vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
    vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
    vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
    vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
    vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
    vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
    vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
    vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    
    /* Clamp to image boundary if outside image */
    if(xn[0]<0) { xn[0]=0;if(xn[1]<0) { xn[1]=0;if(xn[2]<0) { xn[2]=0; if(xn[3]<0) { xn[3]=0; }}}}
    if(yn[0]<0) { yn[0]=0;if(yn[1]<0) { yn[1]=0;if(yn[2]<0) { yn[2]=0; if(yn[3]<0) { yn[3]=0; }}}}
    b=Isize[0]-1;
    if(xn[3]>b) { xn[3]=b;if(xn[2]>b) { xn[2]=b;if(xn[1]>b) { xn[1]=b; if(xn[0]>b) { xn[0]=b; }}}}
    b=Isize[1]-1;
    if(yn[3]>b) { yn[3]=b;if(yn[2]>b) { yn[2]=b;if(yn[1]>b) { yn[1]=b; if(yn[0]>b) { yn[0]=b; }}}}
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for (rgb=0; rgb<3; rgb++) {
        Ipixel[rgb]=0;
        for(i=0; i<4; i++) {
            Ipixelx =vector_qx[0]*getcolor_mindex2(xn[0], yn[i], Isize[0], Isize[1], Iin, rgb);
            Ipixelx+=vector_qx[1]*getcolor_mindex2(xn[1], yn[i], Isize[0], Isize[1], Iin, rgb);
            Ipixelx+=vector_qx[2]*getcolor_mindex2(xn[2], yn[i], Isize[0], Isize[1], Iin, rgb);
            Ipixelx+=vector_qx[3]*getcolor_mindex2(xn[3], yn[i], Isize[0], Isize[1], Iin, rgb);
            Ipixel[rgb]+= vector_qy[i]*Ipixelx;
        }
    }
}

void interpolate_2d_linear_color_black(double *Ipixel, double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1;
    double perc[4]={0, 0, 0, 0};
    double xCom, yCom, xComi, yComi;
    
    /* Loop variable  */
    int rgb;
    
    /*  Rounded location  */
    double fTlocalx, fTlocaly;
    
    Ipixel[0]=0; Ipixel[1]=0; Ipixel[2]=0;
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    xBas1=xBas0+1; yBas1=yBas0+1;
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx; yCom=Tlocaly-fTlocaly;
    xComi=(1-xCom); yComi=(1-yCom);
    perc[0]=xComi * yComi; perc[1]=xComi * yCom; perc[2]=xCom * yComi; perc[3]=xCom * yCom;
    
    if((xBas0>=0)&&(xBas0<Isize[0])) {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            for (rgb=0; rgb<3; rgb++) {
                Ipixel[rgb]+=getcolor_mindex2(xBas0, yBas0, Isize[0], Isize[1], Iin, rgb)*perc[0];
            }
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            for (rgb=0; rgb<3; rgb++) {
                Ipixel[rgb]+=getcolor_mindex2(xBas0, yBas1, Isize[0], Isize[1], Iin, rgb)*perc[1];
            }
        }
    }
    if((xBas1>=0)&&(xBas1<Isize[0]))  {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            for (rgb=0; rgb<3; rgb++) {
                Ipixel[rgb]+=getcolor_mindex2(xBas1, yBas0, Isize[0], Isize[1], Iin, rgb)*perc[2];
            }
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            for (rgb=0; rgb<3; rgb++) {
                Ipixel[rgb]+=getcolor_mindex2(xBas1, yBas1, Isize[0], Isize[1], Iin, rgb)*perc[3];
            }
        }
    }
}

void interpolate_2d_linear_color(double *Ipixel, double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1;
    double perc[4]={0, 0, 0, 0};
    double xCom, yCom, xComi, yComi;
    double color[4]={0, 0, 0, 0};
    
    /* Loop variable  */
    int rgb;
    
    /*  Rounded location  */
    double fTlocalx, fTlocaly;
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    xBas1=xBas0+1; yBas1=yBas0+1;
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx; yCom=Tlocaly-fTlocaly;
    xComi=(1-xCom); yComi=(1-yCom);
    perc[0]=xComi * yComi;
    perc[1]=xComi * yCom;
    perc[2]=xCom * yComi;
    perc[3]=xCom * yCom;
    
    if(xBas0<0) { xBas0=0; if(xBas1<0) { xBas1=0; }}
    if(yBas0<0) { yBas0=0; if(yBas1<0) { yBas1=0; }}
    if(xBas1>(Isize[0]-1)) { xBas1=Isize[0]-1; if(xBas0>(Isize[0]-1)) { xBas0=Isize[0]-1; }}
    if(yBas1>(Isize[1]-1)) { yBas1=Isize[1]-1; if(yBas0>(Isize[1]-1)) { yBas0=Isize[1]-1; }}
    
    for (rgb=0; rgb<3; rgb++) {
        color[0]=getcolor_mindex2(xBas0, yBas0, Isize[0], Isize[1], Iin, rgb);
        color[1]=getcolor_mindex2(xBas0, yBas1, Isize[0], Isize[1], Iin, rgb);
        color[2]=getcolor_mindex2(xBas1, yBas0, Isize[0], Isize[1], Iin, rgb);
        color[3]=getcolor_mindex2(xBas1, yBas1, Isize[0], Isize[1], Iin, rgb);
        Ipixel[rgb]=color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
    }
}

double interpolate_3d_linear_black(double Tlocalx, double Tlocaly, double Tlocalz, int *Isize, double *Iin) {
    double Iout;
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1, zBas0, zBas1;
    double perc[8];
    double xCom, yCom, zCom;
    double xComi, yComi, zComi;
    double color[8]={0, 0, 0, 0, 0, 0, 0, 0};
    double fTlocalx, fTlocaly, fTlocalz;
    
    fTlocalx=floor(Tlocalx); fTlocaly=floor(Tlocaly); fTlocalz=floor(Tlocalz);
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    xBas1=xBas0+1;      yBas1=yBas0+1;      zBas1=zBas0+1;
    
    
    color[0]=0; color[1]=0; color[2]=0; color[3]=0;
    color[4]=0; color[5]=0; color[6]=0; color[7]=0;
    
    if((xBas0>=0)&&(xBas0<Isize[0])) {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[0]=getcolor_mindex3(xBas0, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[1]=getcolor_mindex3(xBas0, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[2]=getcolor_mindex3(xBas0, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[3]=getcolor_mindex3(xBas0, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
    }
    if((xBas1>=0)&&(xBas1<Isize[0]))  {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[4]=getcolor_mindex3(xBas1, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[5]=getcolor_mindex3(xBas1, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[6]=getcolor_mindex3(xBas1, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[7]=getcolor_mindex3(xBas1, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
    }
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx;  yCom=Tlocaly-fTlocaly;   zCom=Tlocalz-fTlocalz;
    
    xComi=(1-xCom); yComi=(1-yCom); zComi=(1-zCom);
    perc[0]=xComi * yComi; perc[1]=perc[0] * zCom; perc[0]=perc[0] * zComi;
    perc[2]=xComi * yCom;  perc[3]=perc[2] * zCom; perc[2]=perc[2] * zComi;
    perc[4]=xCom * yComi;  perc[5]=perc[4] * zCom; perc[4]=perc[4] * zComi;
    perc[6]=xCom * yCom;   perc[7]=perc[6] * zCom; perc[6]=perc[6] * zComi;
    
    /* Set the current pixel value */
    Iout =color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3]+color[4]*perc[4]+color[5]*perc[5]+color[6]*perc[6]+color[7]*perc[7];
    return Iout;
}

double interpolate_3d_linear(double Tlocalx, double Tlocaly, double Tlocalz, int *Isize, double *Iin) {
    double Iout;
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1, zBas0, zBas1;
    double perc[8];
    double xCom, yCom, zCom;
    double xComi, yComi, zComi;
    double color[8]={0, 0, 0, 0, 0, 0, 0, 0};
    double fTlocalx, fTlocaly, fTlocalz;
    
    fTlocalx=floor(Tlocalx); fTlocaly=floor(Tlocaly); fTlocalz=floor(Tlocalz);
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    xBas1=xBas0+1;      yBas1=yBas0+1;      zBas1=zBas0+1;
    
    /* Clamp to boundary */
    if(xBas0<0) {xBas0=0; if(xBas1<0) { xBas1=0; }}
    if(yBas0<0) {yBas0=0; if(yBas1<0) { yBas1=0; }}
    if(zBas0<0) {zBas0=0; if(zBas1<0) { zBas1=0; }}
    if(xBas1>(Isize[0]-1)) { xBas1=Isize[0]-1; if(xBas0>(Isize[0]-1)) { xBas0=Isize[0]-1; }}
    if(yBas1>(Isize[1]-1)) { yBas1=Isize[1]-1; if(yBas0>(Isize[1]-1)) { yBas0=Isize[1]-1; }}
    if(zBas1>(Isize[2]-1)) { zBas1=Isize[2]-1; if(zBas0>(Isize[2]-1)) { zBas0=Isize[2]-1; }}
    
    /*  Get intensities */
    color[0]=getcolor_mindex3(xBas0, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[1]=getcolor_mindex3(xBas0, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    color[2]=getcolor_mindex3(xBas0, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[3]=getcolor_mindex3(xBas0, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    color[4]=getcolor_mindex3(xBas1, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[5]=getcolor_mindex3(xBas1, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    color[6]=getcolor_mindex3(xBas1, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[7]=getcolor_mindex3(xBas1, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx;  yCom=Tlocaly-fTlocaly;   zCom=Tlocalz-fTlocalz;
    
    xComi=(1-xCom); yComi=(1-yCom); zComi=(1-zCom);
    perc[0]=xComi * yComi; perc[1]=perc[0] * zCom; perc[0]=perc[0] * zComi;
    perc[2]=xComi * yCom;  perc[3]=perc[2] * zCom; perc[2]=perc[2] * zComi;
    perc[4]=xCom * yComi;  perc[5]=perc[4] * zCom; perc[4]=perc[4] * zComi;
    perc[6]=xCom * yCom;   perc[7]=perc[6] * zCom; perc[6]=perc[6] * zComi;
    
    /* Set the current pixel value */
    Iout =color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3]+color[4]*perc[4]+color[5]*perc[5]+color[6]*perc[6]+color[7]*perc[7];
    return Iout;
}

double interpolate_3d_cubic_black(double Tlocalx, double Tlocaly, double Tlocalz, int *Isize, double *Iin) {
    /* Floor of coordinate */
    double fTlocalx, fTlocaly, fTlocalz;
    /* Zero neighbor */
    int xBas0, yBas0, zBas0;
    /* The location in between the pixels 0..1 */
    double tx, ty, tz;
    /* Neighbor loccations */
    int xn[4], yn[4], zn[4];
    
    /* The vectors */
    double vector_tx[4], vector_ty[4], vector_tz[4];
    double vector_qx[4], vector_qy[4], vector_qz[4];
    /* Interpolated Intensity; */
    double Ipixelx=0, Ipixelxy=0, Ipixelxyz=0;
    /* Loop variable */
    int i, j;
    
    /* Determine of the zero neighbor */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly); fTlocalz = floor(Tlocalz);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly; tz=Tlocalz-fTlocalz;
    
    /* Determine the t vectors */
    vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
    vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
    vector_tz[0]= 0.5; vector_tz[1]= 0.5*tz; vector_tz[2]= 0.5*pow2(tz); vector_tz[3]= 0.5*pow3(tz);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
    vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
    vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
    vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
    vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
    vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
    vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
    vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
    vector_qz[0]= -1.0*vector_tz[1]+2.0*vector_tz[2]-1.0*vector_tz[3];
    vector_qz[1]= 2.0*vector_tz[0]-5.0*vector_tz[2]+3.0*vector_tz[3];
    vector_qz[2]= 1.0*vector_tz[1]+4.0*vector_tz[2]-3.0*vector_tz[3];
    vector_qz[3]= -1.0*vector_tz[2]+1.0*vector_tz[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    zn[0]=zBas0-1; zn[1]=zBas0; zn[2]=zBas0+1; zn[3]=zBas0+2;
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for(j=0; j<4; j++) {
        Ipixelxy=0;
        if((zn[j]>=0)&&(zn[j]<Isize[2])) {
            for(i=0; i<4; i++) {
                Ipixelx=0;
                if((yn[i]>=0)&&(yn[i]<Isize[1])) {
                    if((xn[0]>=0)&&(xn[0]<Isize[0])) {
                        Ipixelx+=vector_qx[0]*getcolor_mindex3(xn[0], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                    if((xn[1]>=0)&&(xn[1]<Isize[0])) {
                        Ipixelx+=vector_qx[1]*getcolor_mindex3(xn[1], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                    if((xn[2]>=0)&&(xn[2]<Isize[0])) {
                        Ipixelx+=vector_qx[2]*getcolor_mindex3(xn[2], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                    if((xn[3]>=0)&&(xn[3]<Isize[0])) {
                        Ipixelx+=vector_qx[3]*getcolor_mindex3(xn[3], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                }
                Ipixelxy+= vector_qy[i]*Ipixelx;
            }
            Ipixelxyz+= vector_qz[j]*Ipixelxy;
        }
    }
    return Ipixelxyz;
}

double interpolate_3d_cubic(double Tlocalx, double Tlocaly, double Tlocalz, int *Isize, double *Iin) {
    /* Floor of coordinate */
    double fTlocalx, fTlocaly, fTlocalz;
    /* Zero neighbor */
    int xBas0, yBas0, zBas0;
    /* The location in between the pixels 0..1 */
    double tx, ty, tz;
    /* Neighbor loccations */
    int xn[4], yn[4], zn[4];
    
    /* The vectors */
    double vector_tx[4], vector_ty[4], vector_tz[4];
    double vector_qx[4], vector_qy[4], vector_qz[4];
    /* Interpolated Intensity; */
    double Ipixelx=0, Ipixelxy=0, Ipixelxyz=0;
    /* Temporary value boundary */
    int b;
    /* Loop variable */
    int i, j;
    
    /* Determine of the zero neighbor */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly); fTlocalz = floor(Tlocalz);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly; tz=Tlocalz-fTlocalz;
    
    /* Determine the t vectors */
    vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
    vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
    vector_tz[0]= 0.5; vector_tz[1]= 0.5*tz; vector_tz[2]= 0.5*pow2(tz); vector_tz[3]= 0.5*pow3(tz);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
    vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
    vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
    vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
    vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
    vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
    vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
    vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
    vector_qz[0]= -1.0*vector_tz[1]+2.0*vector_tz[2]-1.0*vector_tz[3];
    vector_qz[1]= 2.0*vector_tz[0]-5.0*vector_tz[2]+3.0*vector_tz[3];
    vector_qz[2]= 1.0*vector_tz[1]+4.0*vector_tz[2]-3.0*vector_tz[3];
    vector_qz[3]= -1.0*vector_tz[2]+1.0*vector_tz[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    zn[0]=zBas0-1; zn[1]=zBas0; zn[2]=zBas0+1; zn[3]=zBas0+2;
    
    /* Clamp to boundary */
    if(xn[0]<0) { xn[0]=0;if(xn[1]<0) { xn[1]=0;if(xn[2]<0) { xn[2]=0; if(xn[3]<0) { xn[3]=0; }}}}
    if(yn[0]<0) { yn[0]=0;if(yn[1]<0) { yn[1]=0;if(yn[2]<0) { yn[2]=0; if(yn[3]<0) { yn[3]=0; }}}}
    if(zn[0]<0) { zn[0]=0;if(zn[1]<0) { zn[1]=0;if(zn[2]<0) { zn[2]=0; if(zn[3]<0) { zn[3]=0; }}}}
    b=Isize[0]-1;
    if(xn[3]>b) { xn[3]=b;if(xn[2]>b) { xn[2]=b;if(xn[1]>b) { xn[1]=b; if(xn[0]>b) { xn[0]=b; }}}}
    b=Isize[1]-1;
    if(yn[3]>b) { yn[3]=b;if(yn[2]>b) { yn[2]=b;if(yn[1]>b) { yn[1]=b; if(yn[0]>b) { yn[0]=b; }}}}
    b=Isize[2]-1;
    if(zn[3]>b) { zn[3]=b;if(zn[2]>b) { zn[2]=b;if(zn[1]>b) { zn[1]=b; if(zn[0]>b) { zn[0]=b; }}}}
    
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for(j=0; j<4; j++) {
        Ipixelxy=0;
        for(i=0; i<4; i++) {
            Ipixelx=0;
            Ipixelx+=vector_qx[0]*getcolor_mindex3(xn[0], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelx+=vector_qx[1]*getcolor_mindex3(xn[1], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelx+=vector_qx[2]*getcolor_mindex3(xn[2], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelx+=vector_qx[3]*getcolor_mindex3(xn[3], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelxy+= vector_qy[i]*Ipixelx;
        }
        Ipixelxyz+= vector_qz[j]*Ipixelxy;
    }
    return Ipixelxyz;
}

float interpolate_3d_float_linear_black(float Tlocalx, float Tlocaly, float Tlocalz, int *Isize, float *Iin) {
    float Iout;
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1, zBas0, zBas1;
    float perc[8];
    float xCom, yCom, zCom;
    float xComi, yComi, zComi;
    float color[8]={0, 0, 0, 0, 0, 0, 0, 0};
    float fTlocalx, fTlocaly, fTlocalz;
    
    fTlocalx=floorfloat(Tlocalx); fTlocaly=floorfloat(Tlocaly); fTlocalz=floorfloat(Tlocalz);
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    xBas1=xBas0+1;      yBas1=yBas0+1;      zBas1=zBas0+1;
    
    
    color[0]=0; color[1]=0; color[2]=0; color[3]=0;
    color[4]=0; color[5]=0; color[6]=0; color[7]=0;
    
    if((xBas0>=0)&&(xBas0<Isize[0])) {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[0]=getcolor_mindex3_float(xBas0, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[1]=getcolor_mindex3_float(xBas0, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[2]=getcolor_mindex3_float(xBas0, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[3]=getcolor_mindex3_float(xBas0, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
    }
    if((xBas1>=0)&&(xBas1<Isize[0]))  {
        if((yBas0>=0)&&(yBas0<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[4]=getcolor_mindex3_float(xBas1, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[5]=getcolor_mindex3_float(xBas1, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
        if((yBas1>=0)&&(yBas1<Isize[1])) {
            if((zBas0>=0)&&(zBas0<Isize[2])) {
                color[6]=getcolor_mindex3_float(xBas1, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
            }
            if((zBas1>=0)&&(zBas1<Isize[2])) {
                color[7]=getcolor_mindex3_float(xBas1, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
            }
        }
    }
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx;  yCom=Tlocaly-fTlocaly;   zCom=Tlocalz-fTlocalz;
    
    xComi=(1-xCom); yComi=(1-yCom); zComi=(1-zCom);
    perc[0]=xComi * yComi; perc[1]=perc[0] * zCom; perc[0]=perc[0] * zComi;
    perc[2]=xComi * yCom;  perc[3]=perc[2] * zCom; perc[2]=perc[2] * zComi;
    perc[4]=xCom * yComi;  perc[5]=perc[4] * zCom; perc[4]=perc[4] * zComi;
    perc[6]=xCom * yCom;   perc[7]=perc[6] * zCom; perc[6]=perc[6] * zComi;
    
    /* Set the current pixel value */
    Iout =color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3]+color[4]*perc[4]+color[5]*perc[5]+color[6]*perc[6]+color[7]*perc[7];
    return Iout;
}

float interpolate_3d_float_linear(float Tlocalx, float Tlocaly, float Tlocalz, int *Isize, float *Iin) {
    float Iout;
    /*  Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1, zBas0, zBas1;
    float perc[8];
    float xCom, yCom, zCom;
    float xComi, yComi, zComi;
    float color[8]={0, 0, 0, 0, 0, 0, 0, 0};
    float fTlocalx, fTlocaly, fTlocalz;
    
    fTlocalx=floorfloat(Tlocalx); fTlocaly=floorfloat(Tlocaly); fTlocalz=floorfloat(Tlocalz);
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    xBas1=xBas0+1;      yBas1=yBas0+1;      zBas1=zBas0+1;
    
    /* Clamp to boundary */
    if(xBas0<0) {xBas0=0; if(xBas1<0) { xBas1=0; }}
    if(yBas0<0) {yBas0=0; if(yBas1<0) { yBas1=0; }}
    if(zBas0<0) {zBas0=0; if(zBas1<0) { zBas1=0; }}
    if(xBas1>(Isize[0]-1)) { xBas1=Isize[0]-1; if(xBas0>(Isize[0]-1)) { xBas0=Isize[0]-1; }}
    if(yBas1>(Isize[1]-1)) { yBas1=Isize[1]-1; if(yBas0>(Isize[1]-1)) { yBas0=Isize[1]-1; }}
    if(zBas1>(Isize[2]-1)) { zBas1=Isize[2]-1; if(zBas0>(Isize[2]-1)) { zBas0=Isize[2]-1; }}
    
    /*  Get intensities */
    color[0]=getcolor_mindex3_float(xBas0, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[1]=getcolor_mindex3_float(xBas0, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    color[2]=getcolor_mindex3_float(xBas0, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[3]=getcolor_mindex3_float(xBas0, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    color[4]=getcolor_mindex3_float(xBas1, yBas0, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[5]=getcolor_mindex3_float(xBas1, yBas0, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    color[6]=getcolor_mindex3_float(xBas1, yBas1, zBas0, Isize[0], Isize[1], Isize[2], Iin);
    color[7]=getcolor_mindex3_float(xBas1, yBas1, zBas1, Isize[0], Isize[1], Isize[2], Iin);
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx;  yCom=Tlocaly-fTlocaly;   zCom=Tlocalz-fTlocalz;
    
    xComi=(1-xCom); yComi=(1-yCom); zComi=(1-zCom);
    perc[0]=xComi * yComi; perc[1]=perc[0] * zCom; perc[0]=perc[0] * zComi;
    perc[2]=xComi * yCom;  perc[3]=perc[2] * zCom; perc[2]=perc[2] * zComi;
    perc[4]=xCom * yComi;  perc[5]=perc[4] * zCom; perc[4]=perc[4] * zComi;
    perc[6]=xCom * yCom;   perc[7]=perc[6] * zCom; perc[6]=perc[6] * zComi;
    
    /* Set the current pixel value */
    Iout =color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3]+color[4]*perc[4]+color[5]*perc[5]+color[6]*perc[6]+color[7]*perc[7];
    return Iout;
}

float interpolate_3d_float_cubic_black(float Tlocalx, float Tlocaly, float Tlocalz, int *Isize, float *Iin) {
    /* Floor of coordinate */
    float fTlocalx, fTlocaly, fTlocalz;
    /* Zero neighbor */
    int xBas0, yBas0, zBas0;
    /* The location in between the pixels 0..1 */
    float tx, ty, tz;
    /* Neighbor loccations */
    int xn[4], yn[4], zn[4];
    
    /* The vectors */
    float vector_tx[4], vector_ty[4], vector_tz[4];
    float vector_qx[4], vector_qy[4], vector_qz[4];
    /* Interpolated Intensity; */
    float Ipixelx=0, Ipixelxy=0, Ipixelxyz=0;
    /* Loop variable */
    int i, j;
    /* constant 0.5; */
    const float con=0.5;
    
    /* Determine of the zero neighbor */
    fTlocalx=floorfloat(Tlocalx); fTlocaly=floorfloat(Tlocaly); fTlocalz=floorfloat(Tlocalz);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly; tz=Tlocalz-fTlocalz;
    
    /* Determine the t vectors */
    vector_tx[0]= con; vector_tx[1]= con*tx; vector_tx[2]= con*pow2_float(tx); vector_tx[3]= con*pow3_float(tx);
    vector_ty[0]= con; vector_ty[1]= con*ty; vector_ty[2]= con*pow2_float(ty); vector_ty[3]= con*pow3_float(ty);
    vector_tz[0]= con; vector_tz[1]= con*tz; vector_tz[2]= con*pow2_float(tz); vector_tz[3]= con*pow3_float(tz);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= (float)-1.0*vector_tx[1]+(float)2.0*vector_tx[2]-(float)1.0*vector_tx[3];
    vector_qx[1]= (float)2.0*vector_tx[0]-(float)5.0*vector_tx[2]+(float)3.0*vector_tx[3];
    vector_qx[2]= (float)1.0*vector_tx[1]+(float)4.0*vector_tx[2]-(float)3.0*vector_tx[3];
    vector_qx[3]= (float)-1.0*vector_tx[2]+(float)1.0*vector_tx[3];
    vector_qy[0]= -(float)1.0*vector_ty[1]+(float)2.0*vector_ty[2]-(float)1.0*vector_ty[3];
    vector_qy[1]= (float)2.0*vector_ty[0]-(float)5.0*vector_ty[2]+(float)3.0*vector_ty[3];
    vector_qy[2]= (float)1.0*vector_ty[1]+(float)4.0*vector_ty[2]-(float)3.0*vector_ty[3];
    vector_qy[3]= -(float)1.0*vector_ty[2]+(float)1.0*vector_ty[3];
    vector_qz[0]= -(float)1.0*vector_tz[1]+(float)2.0*vector_tz[2]-(float)1.0*vector_tz[3];
    vector_qz[1]= (float)2.0*vector_tz[0]-(float)5.0*vector_tz[2]+(float)3.0*vector_tz[3];
    vector_qz[2]= (float)1.0*vector_tz[1]+(float)4.0*vector_tz[2]-(float)3.0*vector_tz[3];
    vector_qz[3]= -(float)1.0*vector_tz[2]+(float)1.0*vector_tz[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    zn[0]=zBas0-1; zn[1]=zBas0; zn[2]=zBas0+1; zn[3]=zBas0+2;
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for(j=0; j<4; j++) {
        Ipixelxy=0;
        if((zn[j]>=0)&&(zn[j]<Isize[2])) {
            for(i=0; i<4; i++) {
                Ipixelx=0;
                if((yn[i]>=0)&&(yn[i]<Isize[1])) {
                    if((xn[0]>=0)&&(xn[0]<Isize[0])) {
                        Ipixelx+=vector_qx[0]*getcolor_mindex3_float(xn[0], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                    if((xn[1]>=0)&&(xn[1]<Isize[0])) {
                        Ipixelx+=vector_qx[1]*getcolor_mindex3_float(xn[1], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                    if((xn[2]>=0)&&(xn[2]<Isize[0])) {
                        Ipixelx+=vector_qx[2]*getcolor_mindex3_float(xn[2], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                    if((xn[3]>=0)&&(xn[3]<Isize[0])) {
                        Ipixelx+=vector_qx[3]*getcolor_mindex3_float(xn[3], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
                    }
                }
                Ipixelxy += vector_qy[i]*Ipixelx;
            }
            Ipixelxyz += vector_qz[i]*Ipixelxy;
        }
    }
    return Ipixelxyz;
}
float interpolate_3d_float_cubic(float Tlocalx, float Tlocaly, float Tlocalz, int *Isize, float *Iin) {
    /* Floor of coordinate */
    float fTlocalx, fTlocaly, fTlocalz;
    /* Zero neighbor */
    int xBas0, yBas0, zBas0;
    /* The location in between the pixels 0..1 */
    float tx, ty, tz;
    /* Neighbor loccations */
    int xn[4], yn[4], zn[4];
    
    /* The vectors */
    float vector_tx[4], vector_ty[4], vector_tz[4];
    float vector_qx[4], vector_qy[4], vector_qz[4];
    
    /* Interpolated Intensity; */
    float Ipixelx=0, Ipixelxy=0, Ipixelxyz=0;
    /* Temporary value boundary */
    int b;
    /* Loop variable */
    int i, j;
    /* const 0.5; */
    const float con=0.5;
    
    /* Determine of the zero neighbor */
    fTlocalx = floorfloat(Tlocalx); fTlocaly = floorfloat(Tlocaly); fTlocalz = floorfloat(Tlocalz);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly; tz=Tlocalz-fTlocalz;
    
    /* Determine the t vectors */
    vector_tx[0]= con; vector_tx[1]= con*tx; vector_tx[2]= con*pow2_float(tx); vector_tx[3]= con*pow3_float(tx);
    vector_ty[0]= con; vector_ty[1]= con*ty; vector_ty[2]= con*pow2_float(ty); vector_ty[3]= con*pow3_float(ty);
    vector_tz[0]= con; vector_tz[1]= con*tz; vector_tz[2]= con*pow2_float(tz); vector_tz[3]= con*pow3_float(tz);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= (float)-1.0*vector_tx[1]+(float)2.0*vector_tx[2]-(float)1.0*vector_tx[3];
    vector_qx[1]= (float)2.0*vector_tx[0]-(float)5.0*vector_tx[2]+(float)3.0*vector_tx[3];
    vector_qx[2]= (float)1.0*vector_tx[1]+(float)4.0*vector_tx[2]-(float)3.0*vector_tx[3];
    vector_qx[3]= (float)-1.0*vector_tx[2]+(float)1.0*vector_tx[3];
    vector_qy[0]= -(float)1.0*vector_ty[1]+(float)2.0*vector_ty[2]-(float)1.0*vector_ty[3];
    vector_qy[1]= (float)2.0*vector_ty[0]-(float)5.0*vector_ty[2]+(float)3.0*vector_ty[3];
    vector_qy[2]= (float)1.0*vector_ty[1]+(float)4.0*vector_ty[2]-(float)3.0*vector_ty[3];
    vector_qy[3]= -(float)1.0*vector_ty[2]+(float)1.0*vector_ty[3];
    vector_qz[0]= -(float)1.0*vector_tz[1]+(float)2.0*vector_tz[2]-(float)1.0*vector_tz[3];
    vector_qz[1]= (float)2.0*vector_tz[0]-(float)5.0*vector_tz[2]+(float)3.0*vector_tz[3];
    vector_qz[2]= (float)1.0*vector_tz[1]+(float)4.0*vector_tz[2]-(float)3.0*vector_tz[3];
    vector_qz[3]= -(float)1.0*vector_tz[2]+(float)1.0*vector_tz[3];
    
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    zn[0]=zBas0-1; zn[1]=zBas0; zn[2]=zBas0+1; zn[3]=zBas0+2;
    
    /* Clamp to boundary */
    
    /* Clamp to boundary */
    if(xn[0]<0) { xn[0]=0;if(xn[1]<0) { xn[1]=0;if(xn[2]<0) { xn[2]=0; if(xn[3]<0) { xn[3]=0; }}}}
    if(yn[0]<0) { yn[0]=0;if(yn[1]<0) { yn[1]=0;if(yn[2]<0) { yn[2]=0; if(yn[3]<0) { yn[3]=0; }}}}
    if(zn[0]<0) { zn[0]=0;if(zn[1]<0) { zn[1]=0;if(zn[2]<0) { zn[2]=0; if(zn[3]<0) { zn[3]=0; }}}}
    b=Isize[0]-1;
    if(xn[3]>b) { xn[3]=b;if(xn[2]>b) { xn[2]=b;if(xn[1]>b) { xn[1]=b; if(xn[0]>b) { xn[0]=b; }}}}
    b=Isize[1]-1;
    if(yn[3]>b) { yn[3]=b;if(yn[2]>b) { yn[2]=b;if(yn[1]>b) { yn[1]=b; if(yn[0]>b) { yn[0]=b; }}}}
    b=Isize[2]-1;
    if(zn[3]>b) { zn[3]=b;if(zn[2]>b) { zn[2]=b;if(zn[1]>b) { zn[1]=b; if(zn[0]>b) { zn[0]=b; }}}}
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for(j=0; j<4; j++) {
        Ipixelxy=0;
        for(i=0; i<4; i++) {
            Ipixelx=0;
            Ipixelx+=vector_qx[0]*getcolor_mindex3_float(xn[0], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelx+=vector_qx[1]*getcolor_mindex3_float(xn[1], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelx+=vector_qx[2]*getcolor_mindex3_float(xn[2], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelx+=vector_qx[3]*getcolor_mindex3_float(xn[3], yn[i], zn[j], Isize[0], Isize[1], Isize[2], Iin);
            Ipixelxy+= vector_qy[i]*Ipixelx;
        }
        Ipixelxyz+= vector_qz[i]*Ipixelxy;
    }
    return Ipixelxyz;
}

double interpolate_2d_double_gray(double Tlocalx, double Tlocaly, int *Isize, double *Iin, int cubic, int black) {
    double Ipixel;
    if(cubic) {
        if(black) { Ipixel=interpolate_2d_cubic_gray_black(Tlocalx, Tlocaly, Isize, Iin); }
        else { Ipixel=interpolate_2d_cubic_gray(Tlocalx, Tlocaly, Isize, Iin); }
    }
    else {
        if(black) { Ipixel=interpolate_2d_linear_gray_black(Tlocalx, Tlocaly, Isize, Iin); }
        else { Ipixel=interpolate_2d_linear_gray(Tlocalx, Tlocaly, Isize, Iin); }
    }
    return Ipixel;
}
    
void interpolate_2d_double_color(double *Ipixel, double Tlocalx, double Tlocaly, int *Isize, double *Iin, int cubic, int black) {
    if(cubic) {
        if(black) { interpolate_2d_cubic_color_black(Ipixel, Tlocalx, Tlocaly, Isize, Iin);}
        else { interpolate_2d_cubic_color(Ipixel, Tlocalx, Tlocaly, Isize, Iin);}
    }
    else {
        if(black) { interpolate_2d_linear_color_black(Ipixel, Tlocalx, Tlocaly, Isize, Iin);}
        else { interpolate_2d_linear_color(Ipixel, Tlocalx, Tlocaly, Isize, Iin);}		
    }
}

double interpolate_3d_double_gray(double Tlocalx, double Tlocaly, double Tlocalz, int *Isize, double *Iin, int cubic, int black){
    double Ipixel;
    if(cubic) {
        if(black) { Ipixel=interpolate_3d_cubic_black(Tlocalx, Tlocaly, Tlocalz, Isize, Iin);}
        else { Ipixel=interpolate_3d_cubic(Tlocalx, Tlocaly, Tlocalz, Isize, Iin);}
    }
    else {
        if(black) { Ipixel=interpolate_3d_linear_black(Tlocalx, Tlocaly, Tlocalz, Isize, Iin);}
        else { Ipixel=interpolate_3d_linear(Tlocalx, Tlocaly, Tlocalz, Isize, Iin);}
    }
    return Ipixel;
}
float interpolate_3d_float_gray(float Tlocalx, float Tlocaly, float Tlocalz, int *Isize, float  *Iin, int cubic, int black){
    float Ipixel;
    if(cubic) {
        if(black) { Ipixel=interpolate_3d_float_cubic_black(Tlocalx, Tlocaly, Tlocalz, Isize, Iin);}
        else {Ipixel= interpolate_3d_float_cubic(Tlocalx, Tlocaly, Tlocalz, Isize, Iin);}
    }
    else {
        if(black) { Ipixel=interpolate_3d_float_linear_black(Tlocalx, Tlocaly, Tlocalz, Isize, Iin);}
        else { Ipixel=interpolate_3d_float_linear(Tlocalx, Tlocaly, Tlocalz, Isize, Iin); }
    }
    return Ipixel;
}










