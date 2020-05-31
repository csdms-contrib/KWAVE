/*
Copyright (C) 2020 Luke McGuire
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


//#include<omp.h>
#include<math.h>
#include<stdio.h>
#include <time.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

#define PI 3.141592653589793
#define sqrt2 1.414213562373
#define oneoversqrt2 0.707106781186


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(double *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	m += NR_END;
	m -= nrl;

	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	m += NR_END;
	m -= nrl;

	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ****f4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nwl, long nwh)
/* allocate a double 4tensor with range t[nrl..nrh][ncl..nch][ndl..ndh][nwl..nwh] */
{
	long i,j,k, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,nwid=nwh-nwl+1 ;
	double ****t;

	/* allocate pointers to pointers to pointers to rows */
	t=(double ****) malloc((size_t)((nrow+NR_END)*sizeof(double***)));
	if (!t) nrerror("allocation failure 1 in f4tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to pointers to rows and set pointers to them */
	t[nrl]=(double ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));
	if (!t[nrl]) nrerror("allocation failure 2 in f4tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl][ncl]=(double **) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double*)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f4tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl][ndl]=(double *) malloc((size_t)((nrow*ncol*ndep*nwid+NR_END)*sizeof(double)));
	if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in f4tensor()");
	t[nrl][ncl][ndl] += NR_END;
	t[nrl][ncl][ndl] -= nwl;

    for(i=nrl;i<=nrh;i++)
	{
		if (i > nrl)
		{
			t[i] = t[i-1] + ncol ;
		    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		    t[i][ncl][ndl] = t[i-1][ncl][ndl] + ncol*ndep*nwid ;
		}
		for(j=ncl;j<=nch;j++)
		{
			if (j > ncl)
			{
				t[i][j]=t[i][j-1] + ndep ;
				t[i][j][ndl] = t[i][j-1][ndl] + ndep*nwid ;
			}

			for(k=ndl;k<=ndh;k++)
			{
				if (k > ndl) t[i][j][k] = t[i][j][k-1] + nwid ;
			}
		}
	}

	/* return pointer to pointer to array of pointers to rows */
	return t;
}

void free_f4tensor(double ****t, long nrl, long nrh, long ncl, long nch, 
	long ndl, long ndh, long nwl, long nwh)
/* free a double f4tensor allocated by f4tensor() */
{
	free((FREE_ARG) (t[nrl][ncl][ndl]+nwl-NR_END));
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

// Returns random number drawn from uniform distribtion on [0 1]
double ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double gasdev(idum)
int *idum;
{
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;
        double ran3();

        if  (iset == 0) {
                do {
                        v1=2.0*ran3(idum)-1.0;
                        v2=2.0*ran3(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

void indexx(n,arr,indx)
float arr[];
int indx[],n;
{
        unsigned long i,indxt,ir=n,itemp,j,k,l=1;
        int jstack=0,*istack;
        float a;

        istack=ivector(1,NSTACK);
        for (j=1;j<=n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=1;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[l]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP

void initialcondition(double **TOPO,double ***VEL,double ***U,double ***VELOLD,double ***UOLD,double ***UOLD2,double **DEPTHINIT,double **TOPOINIT,double **MAXVEL,double **MAXDEPTH,int **SOLID,double **KS,double **THETA0,double **THETAS,double **HF,double **VINF,double **CHANNEL,double **MANNING0,double **CV,double **QOLD,int nx,int ny,double hfilm)
{   int i,j;
    
    FILE *fr0,*fr1,*fr2,*fr3,*fr4,*fr5,*fr6,*fr7,*fr8,*fr9,*fr10;
    
    fr0=fopen("./topoin","r");
    fr1=fopen("./depthin","r");
    fr2=fopen("./solidin","r");
    fr3=fopen("./ksin","r");
    fr4=fopen("./theta0in","r");
    fr5=fopen("./thetasin","r");
    fr6=fopen("./hfin","r");
    fr7=fopen("./vinfin","r");
    fr8=fopen("./channelin","r");
    fr9=fopen("./manningin","r");
    fr10=fopen("./vegcoverin","r");

    for (i=1;i<=nx;i++)
        for (j=1;j<=ny;j++)
        {
            fscanf(fr0,"%lf",&TOPO[i][j]);     // TOPOGRAPHY
            fscanf(fr1,"%lf",&U[1][i][j]);     // H
            fscanf(fr2,"%d",&SOLID[i][j]);     // MASK FOR SOLID INTERIOR B.C.
            fscanf(fr3,"%lf",&KS[i][j]);       // SATURATED HYDRAULIC CONDUCTIVITY
            fscanf(fr4,"%lf",&THETA0[i][j]);   // INITIAL SOIL MOISTURE
            fscanf(fr5,"%lf",&THETAS[i][j]);   // SOIL MOISTURE AT SATURATION
            fscanf(fr6,"%lf",&HF[i][j]);
            fscanf(fr7,"%lf",&VINF[i][j]);     // INITIAL VOLUME WATER INFILTRATED  
            fscanf(fr8,"%lf",&CHANNEL[i][j]);     // CHANNEL MASK
            fscanf(fr9,"%lf",&MANNING0[i][j]);        // MANNING COEFFICIENT
            fscanf(fr10,"%lf",&CV[i][j]);        // FRACTION VEGETATION COVER

            TOPOINIT[i][j]=TOPO[i][j];
            DEPTHINIT[i][j]=U[1][i][j];
			U[2][i][j]=0;
			U[3][i][j]=0;
			VEL[1][i][j]=0;
			VEL[2][i][j]=0;
            
            UOLD[1][i][j]=U[1][i][j];
            UOLD[2][i][j]=U[2][i][j];
            UOLD[3][i][j]=U[3][i][j];
            UOLD2[1][i][j]=U[1][i][j];
            UOLD2[2][i][j]=U[2][i][j];
            UOLD2[3][i][j]=U[3][i][j];
            VELOLD[1][i][j]=VEL[1][i][j];
            VELOLD[2][i][j]=VEL[2][i][j];
            QOLD[i][j]=0;
            MAXVEL[i][j]=0;
            MAXDEPTH[i][j]=0;
        }

    fclose(fr0);
    fclose(fr1);
    fclose(fr2);
    fclose(fr3);
    fclose(fr4);
    fclose(fr5);
    fclose(fr6);
    fclose(fr7);
    fclose(fr8);
    fclose(fr9);
    fclose(fr10);
}

int main()
{   FILE *fp0,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10;
    FILE *fpm0,*fpm1,*fr7,*fr11,*fr12;
    int i,j,k,nx,ny,cflstat,rind,stagecnt,rnum,mintopoi,mintopoj;
    double **TOPO,**SX,**INITDEPTH,**DEPTHINIT,*INPUT,**TOPOINIT,***UOLD2,**SF,**SLOPE,**SLOPEX,**SLOPEY,**STAGE;
    double ***U,***VEL,***FLUXX,***FLUXY,**SY,***UOLD,***VELOLD,**MAXVEL,**R,**CV,*RAIN,*RAIN1,*RAIN2,**SGNX,**SGNY,**MANNING,**MANNING0;
    double **KS,**HF,**VINF,**THETA0,**THETAS,**INFL,**ZF,**MAXDEPTH,**CHANNEL,**Q,**QOLD;
    double Sf,cstable,oneoverdx,oneoverdx2,maxstep,q,epsilon,h0,a1,a2,d,d84,flowvel,rint,stageint,stageinterval;
    double dischargeoutx1,dischargeouty1,dischargeoutx2,dischargeouty2,dischargeoutx3,dischargeouty3,dischargeoutx4,dischargeouty4;
    double dt,dx,t,tend,printinterval,printint,velmax,roughnessd,timestep,g,courant,hfilm,minstep,maxdepth1,maxdepth2,maxdepth3,maxdepth4,maxuh,maxvh,depthmax;
    double R1,pi,Si,Ki,gi,CS,Di,evap,mintopo;
    int **SOLID,*BNDRY;
    
    // Set number of threads to use if running in parallel
    //omp_set_num_threads(4);
    
    fp0=fopen("./cflstatus","w");
    
    fpm0=fopen("./depthmovie","w");
    fpm1=fopen("./velocitymovie","w");
    
    fr7=fopen("./input","r");
    fr11=fopen("./rain1","r");
    fr12=fopen("./rain2","r");
    
    // Input Parameters
    INPUT=vector(1,13);
    for (i=1;i<=13;i++)
        fscanf(fr7,"%lf",&INPUT[i]);
    
    // Grid and Time Information
    nx=INPUT[1]; 
    ny=INPUT[2];                        // [m]
    dx=INPUT[3];                        // [m]
    tend=INPUT[4];                      // Total Simulation Time [s]
    
    // Parameters for Hydraulic Roughness
	epsilon=INPUT[5];                    // Empirical Exponent in Depth-Dependent Manning 
	h0=INPUT[6];                         // Empirical Parameter in Depth-Dependent Manning 
    d84=INPUT[9];                       // d84 of channel sediment for use in VPE equation [m]
    
    // Parameters For Rutter Interception Model
    pi=INPUT[10];                            // Free Throughfall Coefficient
    Si=INPUT[11];                            // Canopy capacity [m]
    Ki=INPUT[12];                            // Canopy Drainage Rate Coefficient [m/s]
    gi=INPUT[13];                            // Canopy Drainage Rate Exponent [1/m]
    CS=0;
    Di=0;
    evap=0/(3600*1000);                      // Evaporation rate [m]/[s]
   
    // Rainfall Input Data
    rnum=INPUT[7];                         // Total number of entries in input rainfall time series [-]
    rint=INPUT[8];                         // Time between rainfall entries in input rainfall time series [s]
    RAIN=vector(1,rnum);
    RAIN1=vector(1,rnum);
    RAIN2=vector(1,rnum);
    for (i=1;i<=rnum;i++)
    {
        fscanf(fr11,"%lf",&RAIN1[i]);       // Rain at Gauge 1 [m]/[s]
        fscanf(fr12,"%lf",&RAIN2[i]);       // Rain at Gauge 2 [m]/[s]
    }
    
    // Average rainfall from two rain gages
    for (i=1;i<=rnum;i++)
        RAIN[i]=0.5*(RAIN1[i]+RAIN2[i]);

    oneoverdx=1/dx;
    oneoverdx2=1/(dx*dx);
        
    U=f3tensor(1,3,-1,nx+2,-1,ny+2);                    // Matrix of Conserved Variables
    UOLD=f3tensor(1,3,-1,nx+2,-1,ny+2);                 // Matrix of Conserved Variables
    UOLD2=f3tensor(1,3,-1,nx+2,-1,ny+2);                // Matrix of Conserved Variables
    TOPO=matrix(-1,nx+2,-1,ny+2);                       // Topographic Elevation [m]
    VEL=f3tensor(1,2,-1,nx+2,-1,ny+2);                  // Velocity [m]/[s]
    VELOLD=f3tensor(1,2,-1,nx+2,-1,ny+2);               // Velocity [m]/[s]
    MAXVEL=matrix(1,nx,1,ny);                           // Maximum Velocity [m]/[s]
    MAXDEPTH=matrix(1,nx,1,ny);                         // Maximum Depth [m]/[s]
    SF=matrix(1,nx,1,ny);                               // Friction Slope 
    SLOPE=matrix(1,nx,1,ny);                            // Topographic Slope 
	SLOPEX=matrix(0,nx+1,0,ny+1);                       // Absolute Value of Topographic Slope in X-Direction 
	SGNX=matrix(0,nx+1,0,ny+1);                         // Sign of Topographic Slope in X-Direction 
	SLOPEY=matrix(0,nx+1,0,ny+1);                       // Absolute Value of Topographic Slope in Y-Direction
	SGNY=matrix(0,nx+1,0,ny+1);                         // Sign of Topographic Slope in Y-Direction
	MANNING=matrix(0,nx+1,0,ny+1);                      // Manning Friction Factor
    MANNING0=matrix(0,nx+1,0,ny+1);                      // Manning Friction Factor
    TOPOINIT=matrix(1,nx,1,ny);                         // Initial Topography [m]
    DEPTHINIT=matrix(1,nx,1,ny);                        // Initial Fluid Depth [m]
    SOLID=imatrix(1,nx,1,ny);                           // SOLID[i][j]==1 Defines an Interior B.C.
	R=matrix(1,nx,1,ny);                                // Rainfall Rate [m]/[s]
    INFL=matrix(1,nx,1,ny);                             // Infiltration
    KS=matrix(1,nx,1,ny);                               // Saturated Hydraulic Conductivy
    THETA0=matrix(1,nx,1,ny);                           // Initial Soil Moisture
    THETAS=matrix(1,nx,1,ny);                           // Saturated Zone Soil Moisture
    HF=matrix(1,nx,1,ny);                               // Wetting Front Capillary Pressure Head
    ZF=matrix(1,nx,1,ny);                               // Depth of Wetting Front [m]
    CHANNEL=matrix(1,nx,1,ny);                          // Mask Indicating Channel Locations
    VINF=matrix(1,nx,1,ny);                             // Cumulative Infiltrated Depth [m]
    STAGE=matrix(1,18,1,30000);
    Q=matrix(1,nx,1,ny);                                // Total Runoff Volume From CN-Based Approach
    QOLD=matrix(1,nx,1,ny);                             // Total Runoff Volume From CN-Based Approach
    R=matrix(1,nx,1,ny);                                // Rainfall Intensity [m]/[s]
    CV=matrix(1,nx,1,ny);                               // Fraction Vegetation Cover within Grid Cell [-]
    
    FLUXX=f3tensor(1,3,0,nx+1,0,ny);            // X-direction Flux
    FLUXY=f3tensor(1,3,0,nx,0,ny+1);            // Y-direction Flux
    
    // Boudnary Conditions
    BNDRY=ivector(1,4);             // 0==Transmissive, 1==Solid
    BNDRY[1]=0;                     // B.C. at i==1
    BNDRY[2]=0;                     // B.C. at j==1      
    BNDRY[3]=0;                     // B.C. at i==nx
    BNDRY[4]=0;                     // B.C. at j==ny
    
    // Time Parameters
    t=0;                            // [s]
    dt=0.0001;                      // Initial Time Step for Fluid Equations [s]
    printinterval=fmax(1,tend/5);
    printint=fmax(1,tend/5);
    stageint=5;
    stageinterval=5;
    stagecnt=0;
    cflstat=0;                                  // clfstat==1 indicates cfl error
    cstable=0.2;                                // Adaptive time step will adjust to make courant=cstable
    maxstep=0.1;                                // Maximum Time Step Regardless of Courant Number
    minstep=1e-5;                               // Simulation will abort if time step is too small
    
    // Additional Model Parameters
    g=9.81;                         // [m]/[s^2]
    hfilm=1e-5;                     // Water Depth < hfilm is Treated as Dry
    a1=6.5;                         // Empirical coefficient used to compute manning n (value from Furgeson [2007])
    a2=2.5;                         // Empirical coefficient used to compute manning n (value from Furgeson [2007])   

    // Set Initial Conditions
    initialcondition(TOPO,VEL,U,VELOLD,UOLD,UOLD2,DEPTHINIT,TOPOINIT,MAXVEL,MAXDEPTH,SOLID,KS,THETA0,THETAS,HF,VINF,CHANNEL,MANNING0,CV,QOLD,nx,ny,hfilm);
    
    mintopo=99999;
	// Compute Topographic Slope
	for (i=1;i<=nx;i++)
		for (j=1;j<=ny;j++)
		{            
            if (TOPO[i][j]<mintopo && SOLID[i][j]==0)
            {
                mintopo=TOPO[i][j];
                mintopoi=i;
                mintopoj=j;
            }
			if (SOLID[i][j]==0)
			{
				// Compute Slope
				if (i!=1 && i!=nx) 
				{
					SGNX[i][j]=1;
					SLOPEX[i][j]=fabs(0.5*oneoverdx*(TOPO[i+1][j]-TOPO[i-1][j]));
					if ((TOPO[i+1][j]-TOPO[i-1][j])<0) SGNX[i][j]=-1;
					if (SOLID[i+1][j]==1) 
					{
						SLOPEX[i][j]=fabs(oneoverdx*(TOPO[i][j]-TOPO[i-1][j])); 
						if ((TOPO[i][j]-TOPO[i-1][j])<0) SGNX[i][j]=-1;
					}
					if (SOLID[i-1][j]==1) 
					{
						SLOPEX[i][j]=fabs(oneoverdx*(TOPO[i+1][j]-TOPO[i][j])); 
						if ((TOPO[i+1][j]-TOPO[i][j])<0) SGNX[i][j]=-1;
					}
				}
				if (j!=1 && j!=ny) 
				{
					SGNY[i][j]=1;
					SLOPEY[i][j]=fabs(0.5*oneoverdx*(TOPO[i][j+1]-TOPO[i][j-1]));
					if ((TOPO[i][j+1]-TOPO[i][j-1])<0) SGNY[i][j]=-1;
					if (SOLID[i][j+1]==1)
					{
						SLOPEY[i][j]=fabs(oneoverdx*(TOPO[i][j]-TOPO[i][j-1]));
						if ((TOPO[i][j]-TOPO[i][j-1])<0) SGNY[i][j]=-1;
					}
					if (SOLID[i][j-1]==1)
					{
						SLOPEY[i][j]=fabs(oneoverdx*(TOPO[i][j+1]-TOPO[i][j]));
						if ((TOPO[i][j+1]-TOPO[i][j])<0) SGNY[i][j]=-1;
					}
				}
					
				if (i==1) 
				{
					SGNX[i][j]=1;
					SLOPEX[i][j]=fabs(oneoverdx*(TOPO[i+1][j]-TOPO[i][j]));
					if ((TOPO[i+1][j]-TOPO[i][j])<0) SGNX[i][j]=-1;
				}
				if (i==nx) 
				{
					SGNX[i][j]=1;
					SLOPEX[i][j]=fabs(oneoverdx*(TOPO[i][j]-TOPO[i-1][j]));
					if ((TOPO[i][j]-TOPO[i-1][j])<0) SGNX[i][j]=-1;
				}
				if (j==1) 
				{
					SGNY[i][j]=1;
					SLOPEY[i][j]=fabs(oneoverdx*(TOPO[i][j+1]-TOPO[i][j]));
					if ((TOPO[i][j+1]-TOPO[i][j])<0) SGNY[i][j]=-1;
				}
				if (j==ny) 
				{
					SGNY[i][j]=1;
					SLOPEY[i][j]=fabs(oneoverdx*(TOPO[i][j]-TOPO[i][j-1]));
					if ((TOPO[i][j]-TOPO[i][j-1])<0) SGNY[i][j]=-1;
				}
                
                if (SLOPEY[i][j]<1e-3 && SLOPEX[i][j]<1e-3)
                {
                    printf("!! Error: Slope is Approximately Zero at (i,j)=(%d,%d)\n",i,j);
                    printf("!! TOPO at (i,j)=%lf\n",TOPO[i][j]);
                    printf("!! SLOPEX at (i,j)=%lf\n",SLOPEX[i][j]);
                    printf("!! SLOPEY at (i,j)=%lf\n",SLOPEY[i][j]);
                }
			}
			else 
			{
				SGNX[i][j]=1;
				SLOPEX[i][j]=1e-3;
				SGNY[i][j]=1;
				SLOPEY[i][j]=1e-3;
			}
            SLOPE[i][j]=pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j],0.5);
		}

    printf("Channel Outlet Elevation [m]: %lf\nOutlet Pixel: (%d,%d)\n",TOPO[mintopoi][mintopoj],mintopoi,mintopoj);
    
	// Topography Boundary Conditions
	for (j=0;j<=ny;j++)
	{
		SLOPEX[0][j]=SLOPEX[1][j];
		SGNX[0][j]=SGNX[1][j];
		SLOPEX[nx+1][j]=SLOPEX[nx][j];
		SGNX[nx+1][j]=SGNX[nx][j];
	}
	for (i=0;i<=nx;i++)
	{
		SLOPEY[i][0]=SLOPEY[i][1];
		SGNY[i][0]=SGNY[i][1];
		SLOPEY[i][ny+1]=SLOPEY[i][ny];
        SGNY[i][ny+1]=SGNY[i][ny];
    }
    
    // Flux Boundary Conditions
    // No Water Enters The Domain Through Boundary
    for (j=0;j<=ny+1;j++)
    {
        FLUXX[1][0][j]=0;
        FLUXX[1][nx+1][j]=0;
        FLUXY[1][0][j]=0;
        FLUXY[1][nx+1][j]=0;
    }
    for (i=0;i<nx+1;i++)
    {
        FLUXX[1][i][0]=0;
        FLUXX[1][i][ny+1]=0;
        FLUXY[1][i][0]=0;
        FLUXY[1][i][ny+1]=0;
    }
    
    time_t t0 = time(NULL);
    
    t=0;
    while (t<tend)
    {
        velmax=0;
        depthmax=0;
        
        // Determine Rainfall Rate From Input Data
        rind=fmin(ceil(t/rint+0.01),rnum);
        R1=RAIN[rind];                                  // Rainfall rate [m]/[s]
        
        // Compute Fluxes
        //#pragma omp parallel for private(j,velmax,depthmax)
        for (i=1;i<=nx;i++)
            for (j=1;j<=ny;j++)
            {
                // Rutter Interception Model
                R[i][j]=(1-CV[i][j])*R1+CV[i][j]*(pi*R1+Di);                      // Effective Rainfall Rate [m/s]
                
                // Manning Friciton Factor
				if (U[1][i][j]>h0) MANNING[i][j]=MANNING0[i][j];
				else MANNING[i][j]=MANNING0[i][j]*pow(U[1][i][j]/h0,epsilon);
                
                // Flux In X-Direction
                if (U[1][i][j]>hfilm)
                {
                    // Compute Velocity using Manning 
                    VEL[1][i][j]=-SGNX[i][j]*pow(U[1][i][j],0.667)*pow(SLOPEX[i][j],0.5)/MANNING[i][j];
                    
                    // Compute Velocity using VPE
                    //if (CHANNEL[i][j]==1) VEL[1][i][j]=-SGNX[i][j]*pow(g*U[1][i][j]*SLOPEX[i][j],0.5)*(a1*U[1][i][j]/d84)/pow((pow(U[1][i][j]/d84,1.6667)+pow(a1/a2,2)),0.5);
                    //else VEL[1][i][j]=-SGNX[i][j]*pow(U[1][i][j],0.667)*pow(SLOPEX[i][j],0.5)/MANNING[i][j];
                }
                else
                {
                    VEL[1][i][j]=0;
                }
                FLUXX[1][i][j]=U[1][i][j]*VEL[1][i][j];
                
                // Flux In Y-Direction
                if (U[1][i][j]>hfilm)
                {
                    // Compute Velocity using Manning 
                    VEL[2][i][j]=-SGNY[i][j]*pow(U[1][i][j],0.667)*pow(SLOPEY[i][j],0.5)/MANNING[i][j];
                    
                    // Compute Velocity using VPE
                    //if (CHANNEL[i][j]==1) VEL[2][i][j]=-SGNY[i][j]*pow(g*U[1][i][j]*SLOPEY[i][j],0.5)*(a1*U[1][i][j]/d84)/pow((pow(U[1][i][j]/d84,1.6667)+pow(a1/a2,2)),0.5);
                    //else VEL[2][i][j]=-SGNY[i][j]*pow(U[1][i][j],0.667)*pow(SLOPEY[i][j],0.5)/MANNING[i][j];
                }
                else
                {
                    VEL[2][i][j]=0;
                }
                FLUXY[1][i][j]=U[1][i][j]*VEL[2][i][j];
                
                // Interior Boundary Conditions
                if (SOLID[i][j]==1)
                {
                    FLUXX[1][i][j]=0;
                    FLUXY[1][i][j]=0;
                }
                
                if (VEL[1][i][j]>velmax) velmax=VEL[1][i][j];
                if (VEL[2][i][j]>velmax) velmax=VEL[2][i][j];
                if (U[1][i][j]>depthmax) depthmax=U[1][i][j];
            }
        
        // Adaptive Time Step
        dt=cstable*dx/(velmax);
        if (dt>maxstep) dt=maxstep;
        if (t+dt>tend) dt=tend-t;
        courant=velmax*dt*oneoverdx;
        if (courant>1)  {printf("!!! CFL Error at t=%lf\n",t); t=2*tend;}
        
        // Rutter Interception Model
        Di=0;
        if (CS>Si) Di=Ki*exp(gi*(CS-Si));
        CS=CS+dt*(1-pi)*R1-dt*Di-fmin(CS,dt*evap);        // Canopy Storage
        
        // Time Update
        //#pragma omp parallel for private(j,flowvel)
        for (i=1;i<=nx;i++)
            for (j=1;j<=ny;j++)
            {
                
                if (SOLID[i][j]==0)
                {
                    // Time Update
                    U[1][i][j] = U[1][i][j]
                            - dt*oneoverdx*(fabs(FLUXX[1][i][j])+fabs(FLUXY[1][i][j]))                  // Flux of Water Out of Cell
                            + dt*oneoverdx*(fmax(FLUXX[1][i-1][j],0)-fmin(FLUXX[1][i+1][j],0))  // Flux of Water Into Cell From X-Direction
                            + dt*oneoverdx*(fmax(FLUXY[1][i][j-1],0)-fmin(FLUXY[1][i][j+1],0))  // Flux of Water Into Cell From Y-Direction
                            + dt*R[i][j];
                }
                else
                {
                    // Pixel Is Not In Computational Domain
                    U[1][i][j]=0;
                    VEL[1][i][j]=0;
                    VEL[2][i][j]=0;
                }
                
                // Infiltration
                ZF[i][j]=VINF[i][j]/(THETAS[i][j]-THETA0[i][j]);                                         // Depth of wetting front
                INFL[i][j]=fmin(U[1][i][j]-hfilm,dt*KS[i][j]*(ZF[i][j]+HF[i][j]+U[1][i][j])/ZF[i][j]);   // Infiltration Capacity
                                        
                VINF[i][j]+=INFL[i][j];                                                                 // Infiltrated depth
                U[1][i][j]-=INFL[i][j];
                QOLD[i][j]=Q[i][j];
                                
                flowvel=pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5);
                if (flowvel>MAXVEL[i][j]) MAXVEL[i][j]=flowvel;
                if (U[1][i][j]>MAXDEPTH[i][j]) MAXDEPTH[i][j]=U[1][i][j];
                
                if (U[1][i][j]<0)
                {
                    printf("!! Error: Flow Depth < 0 at (i=%d,j=%d)\n",i,j);
                    t=tend+1;
                }
                
                U[2][i][j]=VEL[1][i][j]*U[1][i][j];
                U[3][i][j]=VEL[2][i][j]*U[1][i][j];
                
            }

        t+=dt;
        
        // Print Intermediate Results
        if (t>printint)
        {
            printf("Time: %lf\t dt: %lf\t R (mm/h):%lf\n",t,dt,R[1][1]*3600*1000);
            printint+=printinterval;
            for (j=1;j<=ny;j++)
                for (i=1;i<=nx;i++)
                {
                    if (i!=nx)
                    {
                        fprintf(fpm0,"%10.10lf\t",U[1][i][j]);
                        fprintf(fpm1,"%10.10lf\t",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                    }
                    if (i==nx)
                    {
                        fprintf(fpm0,"%10.10lf\n",U[1][i][j]);
                        fprintf(fpm1,"%10.10lf\n",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                    }
                }
        }
        
        if (t>stageint)
        {
            stageint+=stageinterval;
            stagecnt++;
           
            STAGE[1][stagecnt]=t;
            
            dischargeoutx1=0;
            dischargeouty1=0;
            dischargeoutx2=0;
            dischargeouty2=0;
            dischargeoutx3=0;
            dischargeouty3=0;
            dischargeoutx4=0;
            dischargeouty4=0;
            maxdepth1=0;
            maxdepth2=0;
            maxdepth3=0;
            maxdepth4=0;
            maxuh=0;
            maxvh=0;
            
            // discharge data at x boundaries
            for (j=1;j<=ny;j++)
            {
                dischargeouty1+=fabs(U[1][nx][j]*VEL[2][nx][j]);
                dischargeoutx1+=fabs(U[1][nx][j]*VEL[1][nx][j]);
                dischargeouty2+=fabs(U[1][1][j]*VEL[2][1][j]);
                dischargeoutx2+=fabs(U[1][1][j]*VEL[1][1][j]);
                if (U[1][nx][j]>maxdepth1) maxdepth1=U[1][nx][j];
                if (U[1][1][j]>maxdepth2) maxdepth2=U[1][1][j];
                if (fabs(U[1][nx][j]*VEL[2][nx][j])>maxuh) maxuh=fabs(U[1][nx][j]*VEL[2][nx][j]);
                if (fabs(U[1][nx][j]*VEL[1][nx][j])>maxvh) maxvh=fabs(U[1][nx][j]*VEL[1][nx][j]);
            }
            
            // discharge data at y boundaries
            for (i=1;i<=nx;i++)
            {
                dischargeouty3+=fabs(U[1][i][ny]*VEL[2][i][ny]);
                dischargeoutx3+=fabs(U[1][i][ny]*VEL[1][i][ny]);
                dischargeouty4+=fabs(U[1][i][1]*VEL[2][i][1]);
                dischargeoutx4+=fabs(U[1][i][1]*VEL[1][i][1]);
                if (U[1][i][ny]>maxdepth3) maxdepth3=U[1][i][ny];
                if (U[1][i][1]>maxdepth4) maxdepth4=U[1][i][1];
            }
            
            STAGE[2][stagecnt]=maxuh;
            STAGE[3][stagecnt]=maxvh;
            STAGE[4][stagecnt]=dischargeoutx1;
            STAGE[5][stagecnt]=dischargeouty1;
            STAGE[6][stagecnt]=dischargeoutx2;
            STAGE[7][stagecnt]=dischargeouty2;
            STAGE[8][stagecnt]=dischargeoutx3;
            STAGE[9][stagecnt]=dischargeouty3;
            STAGE[10][stagecnt]=dischargeoutx4;
            STAGE[11][stagecnt]=dischargeouty4;
            STAGE[12][stagecnt]=maxdepth1;
            STAGE[13][stagecnt]=maxdepth2;
            STAGE[14][stagecnt]=maxdepth3;
            STAGE[15][stagecnt]=maxdepth4;
            STAGE[16][stagecnt]=U[1][mintopoi][mintopoj];
            STAGE[17][stagecnt]=U[1][mintopoi][mintopoj]*fabs(VEL[1][mintopoi][mintopoj]);
            STAGE[18][stagecnt]=U[1][mintopoi][mintopoj]*fabs(VEL[2][mintopoi][mintopoj]);
        }

    }
    
    time_t t1 = time(NULL);
    printf ("Elapsed wall clock time: %ld\n", (long) (t1 - t0));
    
    // Print Final Results
    fprintf(fp0,"%d\n",cflstat);
    
    fp1=fopen("./topo","w");
    fp2=fopen("./depth","w");
    fp3=fopen("./u","w");
    fp4=fopen("./v","w");
    fp5=fopen("./vel","w");
    fp6=fopen("./solid","w");
    fp7=fopen("./maxvel","w");
    fp8=fopen("./maxdepth","w");
    fp9=fopen("./infl","w");
    fp10=fopen("./stage","w");
    
    for (j=1;j<=ny;j++)
        for (i=1;i<=nx;i++)
        {
            if (i!=nx)
            {
                fprintf(fp1,"%10.10lf\t",TOPO[i][j]);
                fprintf(fp2,"%10.10lf\t",U[1][i][j]);
                fprintf(fp3,"%10.10lf\t",VEL[1][i][j]);
                fprintf(fp4,"%10.10lf\t",VEL[2][i][j]);
                //fprintf(fp3,"%10.10lf\t",SLOPEX[i][j]);
                //fprintf(fp4,"%10.10lf\t",SLOPEY[i][j]);
                fprintf(fp5,"%10.10lf\t",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                fprintf(fp6,"%d\t",SOLID[i][j]);
                fprintf(fp7,"%10.10lf\t",MAXVEL[i][j]);
                fprintf(fp8,"%10.10lf\t",MAXDEPTH[i][j]);
                fprintf(fp9,"%10.10lf\t",VINF[i][j]);
            }
            if (i==nx)
            {
                fprintf(fp1,"%10.10lf\n",TOPO[i][j]);
                fprintf(fp2,"%10.10lf\n",U[1][i][j]);
                fprintf(fp3,"%10.10lf\n",VEL[1][i][j]);
                fprintf(fp4,"%10.10lf\n",VEL[2][i][j]);
                //fprintf(fp3,"%10.10lf\n",SLOPEX[i][j]);
                //fprintf(fp4,"%10.10lf\n",SLOPEY[i][j]);
                fprintf(fp5,"%10.10lf\n",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                fprintf(fp6,"%d\n",SOLID[i][j]);
                fprintf(fp7,"%10.10lf\n",MAXVEL[i][j]);
                fprintf(fp8,"%10.10lf\n",MAXDEPTH[i][j]);
                fprintf(fp9,"%10.10lf\n",VINF[i][j]);
            }
        }
    
    for (j=1;j<=stagecnt;j++)
        for (i=1;i<=18;i++)
        {
            if (i!=18) fprintf(fp10,"%lf\t",STAGE[i][j]);
            if (i==18) fprintf(fp10,"%lf\n",STAGE[i][j]);
        }
    
    
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    fclose(fp6);
    fclose(fp7);
    fclose(fp8);
    fclose(fp9);
    fclose(fp10);
    
    fclose(fpm0);
    fclose(fpm1);
    
    fclose(fr7);
    
}
