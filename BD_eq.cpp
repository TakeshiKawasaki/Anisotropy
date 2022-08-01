# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <time.h>
using namespace std;
# define Npm 16000
# define NR 20000
# define Nn 50
# define Pm 500
# define pi 3.141592653589793238462643383279
int copy(double *x_update,double *y_update,double *x,double *y,int Np,double x_corr,double y_corr){
  int i;
  for(i=0;i<Np;i++){
    x_update[i]=x[i]-x_corr;
    y_update[i]=y[i]-y_corr;
  }
  return 0;
}

int  set_map(double *x,double *y,int (*map)[Nn],double *a,double L,int Np,double chi,double *theta)
{
  int i,j;
  int n;
  double r2;
  double aij;
  double dx,dy;
  double Ai,Aj,A;
  
  for(i=0;i<Np;i++)
    {
      n=0;
      for(j=0;j<Np;j++)
	if(i!=j)
	  {
	    dx=x[i]-x[j];
	    dy=y[i]-y[j];
	    
	    if(dx>L/2)
	      dx-=L;
	    if(dx<-L/2)
	      dx+=L;
	    
	    if(dy> L/2.)
	      dy-=L;
	    if(dy<-L/2.)
	      dy+=L;

	    r2=dx*dx+dy*dy;
	    Ai=cos(theta[i])*dx+sin(theta[i])*dy;
	    Aj=cos(theta[j])*dx+sin(theta[j])*dy;
	    A=  chi*(Ai*Ai+Aj*Aj)/r2;
	    aij=0.5*(a[i]+a[j]);     
	    
	    if(r2< aij*aij*pow((1.+A),1./3.)*1.2*1.2)
	      {
		// cout << pow((1.+A),1./6.)<< endl;
		n++;
		map[i][n]=j;
	      }
	  }
      map[i][0]=n;
    }
  
  return 0;
}

int sq(double *x,double *y,double *qx,double *qy,double L,int Np){
  int i,j,nx,ny;
  double si=0.,co=0.,c2,s2,temp,sq[1000000],tmp;
  int nmax=150;
  
  for(nx=-nmax;nx<nmax;nx++)
    for(ny=-nmax;ny<nmax;ny++){
      co=0.0;
      si=0.0;      
      for(i=0;i<Np;i++){
	///(0.5*pi*int(sqrt(nx*nx+ny*ny)))/Np/ens
	co+=cos(2.*pi*(nx*x[i]+ny*y[i])/L);
	si+=sin(2.*pi*(nx*x[i]+ny*y[i])/L);
      }        
      c2=co*co;
      s2=si*si;
      sq[nmax+nx+2*nmax*(ny+nmax)]=(c2+s2)/Np;
      qx[nmax+nx+2*nmax*(ny+nmax)]=2.*pi*nx/L;
      qy[nmax+nx+2*nmax*(ny+nmax)]=2.*pi*ny/L;
    }
  
  //sorting
  for (i=0; i<4*nmax*nmax; ++i) {
    for (j=i+1; j<4*nmax*nmax; ++j) {
      if (sq[j] > sq[i]) {
	tmp =  sq[i];
	sq[i] = sq[j];
	sq[j] = tmp;
	
	tmp = qx[i];
	qx[i] = qx[j];
	qx[j] = tmp;
	
	tmp =  qy[i];
	qy[i] = qy[j];
	qy[j] = tmp;
	
      }
    }
  }
  for(i=0;i<10;i++)
    cout<<"sq="<<sq[i] <<" nx="<<qx[i]<<" ny="<<qy[i]<<endl;
  
  // cout <<"sqfinish"<<endl;
  
  return 0;
}



int calc_order(double *x,double *y,double *theta,int Np, double L,double *phi6,double *phi4,double *phi2,double *s,double *avphi6,double *avphi4,double *avphi2,double *avs,double *a,double *gr,double *g6,double *g4,double *g2,double *gs,double chi,double *qx,double *qy,double *cq,double *avphiq){
  int (*map)[Nn]=new int[Npm][Nn];
  double angle,*phi6c,*phi6s,*phi4c,*phi4s,*phi2c,*phi2s,av6s=0.0,av6c=0.0,av4s=0.0,av4c=0.0,av2s=0.0,av2c=0.0,dx,dy,r;
  double t;
  
  phi6c=new double[Npm];
  phi6s=new double[Npm];
  phi4c=new double[Npm];
  phi4s=new double[Npm];
  phi2c=new double[Npm];
  phi2s=new double[Npm];
  double phiqx=0.0;
  double phiqy=0.0;
  set_map(x,y,map,a,L,Np,chi,theta);
  
  int i,j,k;
  // *avphi6=0.0;
  // *avphi4=0.0;
  // *avs=0.0;
  for(i=0;i<Np;i++){
    phi6c[i]=0.0;
    phi6s[i]=0.0;
    phi4c[i]=0.0;
    phi4s[i]=0.0;
    phi2c[i]=0.0;
    phi2s[i]=0.0;
    s[i]=0.0;
    for(j=1;j<=map[i][0];j++){
      dx=x[map[i][j]]-x[i];
      if(dx>0.5*L) dx-=L;
      if(dx<-0.5*L) dx+=L;
      dy=y[map[i][j]]-y[i];
      if(dy>0.5*L) dy-=L;
      if(dy<-0.5*L) dy+=L;
      angle = atan(dx/(dy+0.000000001));

      phi2c[i]+=cos(2.*angle)/(map[i][0]+0.000001);
      phi2s[i]+=sin(2.*angle)/(map[i][0]+0.000001);

      phi4c[i]+=cos(4.*angle)/(map[i][0]+0.000001);
      phi4s[i]+=sin(4.*angle)/(map[i][0]+0.000001);

      phi6c[i]+=cos(6.*angle)/(map[i][0]+0.000001);
      phi6s[i]+=sin(6.*angle)/(map[i][0]+0.000001);
      t=cos(theta[i])*cos(theta[map[i][j]])+sin(theta[i])*sin(theta[map[i][j]]);
      s[i]+=(2.*t*t-1);
    }

    phiqx+=cos(qx[1]*x[i]+qy[1]*y[i]);
    phiqy+=sin(qx[1]*x[i]+qy[1]*y[i]);
     
    s[i]/=(map[i][0]+0.000001);
    *avs+=s[i]/Np;
    phi6[i]=sqrt(phi6c[i]*phi6c[i]+phi6s[i]*phi6s[i]);
    phi4[i]=sqrt(phi4c[i]*phi4c[i]+phi4s[i]*phi4s[i]);
    phi2[i]=sqrt(phi2c[i]*phi2c[i]+phi2s[i]*phi2s[i]);
    av6s+=phi6s[i];
    av6c+=phi6c[i];
    av4s+=phi4s[i];
    av4c+=phi4c[i];
    av2s+=phi2s[i];
    av2c+=phi2c[i];
  }
  *avphiq+=sqrt(phiqx*phiqx+phiqy*phiqy)/Np;
  *avphi6+=sqrt(av6s*av6s+av6c*av6c)/Np;
  *avphi4+=sqrt(av4s*av4s+av4c*av4c)/Np;
  *avphi2+=sqrt(av2s*av2s+av2c*av2c)/Np;

  for(i=0;i<Np;i++){
    for(j=0;j<Np;j++){
      if(i!=j){
	dx=x[j]-x[i];
	if(dx>0.5*L) dx-=L;
	if(dx<-0.5*L) dx+=L;
	dy=y[j]-y[i];
	if(dy>0.5*L) dy-=L;
	if(dy<-0.5*L) dy+=L;
	r=sqrt(dx*dx+dy*dy);
	gr[int(r/0.02)]+=1.0/(2.0*pi*r*0.02*Np*Np/L/L);
	g6[int(r/0.02)]+=(phi6s[i]*phi6s[j]+phi6c[i]*phi6c[j])/(2.0*pi*r*0.02*Np*Np/L/L);
	g4[int(r/0.02)]+=(phi4s[i]*phi4s[j]+phi4c[i]*phi4c[j])/(2.0*pi*r*0.02*Np*Np/L/L);
	g2[int(r/0.02)]+=(phi2s[i]*phi2s[j]+phi2c[i]*phi2c[j])/(2.0*pi*r*0.02*Np*Np/L/L);
	t=cos(theta[i])*cos(theta[j])+sin(theta[i])*sin(theta[j]);
	gs[int(r/0.02)]+=(2.0*t*t-1.0)/(2.0*pi*r*0.02*Np*Np/L/L);

	for(k=1;k<=20;k++)
	  cq[int(r/0.02)]+=cos(qx[k]*dx+qy[k]*dy)/(20.0*2.0*pi*r*0.02*Np*Np/L/L);
      }
    }
  }

  delete [] phi6s;
  delete [] phi6c;
  delete [] phi4s;
  delete [] phi4c;
  delete [] phi2s;
  delete [] phi2c;
  delete [] map;
  return 0;
}

int calc_disp(double *disp_max,double *disp_ave,double *x,double *y,double *x_update,double *y_update,double *x0,double *y0,int Np,int L,double x_corr,double y_corr){
  int i;
  double dx,dy;
  double disp;
  
  *disp_ave = 0.0;

  for(i=0;i<Np;i++){
    dx=x[i]-x_update[i]-x_corr;
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    dy=y[i]-y_update[i]-y_corr;
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    disp = dx*dx+dy*dy;
    if(disp > *disp_max)
      *disp_max =disp;    
  }
  
  for(i=0;i<Np;i++){
    dx=x[i]-x0[i]-x_corr;
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    dy=y[i]-y0[i]-y_corr;
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    disp = dx*dx+dy*dy;
    *disp_ave+=disp/Np;
  }
  
  return 0;
}

int calc_disp_max(double *disp_max,double *x,double *y,double *x_update,double *y_update,int Np,int L){
  int i;
  double dx,dy;
  double disp;
  
  for(i=0;i<Np;i++){
    dx=x[i]-x_update[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    
    dy=y[i]-y_update[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    disp = dx*dx+dy*dy;    
    if(disp > *disp_max)
      *disp_max =disp;    
  }
  
  return 0;
}


int com_correction(double *x,double *y,double *x_corr,double *y_corr,int Np,double L){
  int i;
  double dx,dy;
  static double x0[Npm],y0[Npm];
  static bool IsFirst = true;
  if(IsFirst){
    for(i=0;i<Np;i++){
      x0[i]=x[i];
      y0[i]=y[i];
    }
    IsFirst = false;
  }
  
  for(i=0;i<Np;i++){
    dx=x[i]-x0[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    
    dy=y[i]-y0[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    
    *x_corr+=dx/Np; //center of mass displacement.x
    *y_corr+=dy/Np;
    
    x0[i]=x[i];
    y0[i]=y[i];
  }
  return 0;
}
int p_bound(double* x, double* y, int Np, double L) {
  int k;
  for (k = 0; k < Np; k++) {
    if (x[k] < 0.0) {
      x[k] = x[k] + L;
    }
    if (x[k] >  L) {
      x[k] = x[k] - L;
    }
    if (y[k] < 0.0) {
      y[k] = y[k] + L;
    }    
    if (y[k] >  L) {
      y[k] = y[k] - L;
    }
  }
  return 0;
}

int output(double *x,double *y,double x_corr,double y_corr,double L,double *theta,double *r1,int Np,double t,double dt,double *phi6,double *phi4,double *phi2,double *s,char *argv[],double *gr,double *g6,double *g4,double *gs,double *cq,int count,double avs,int count_th){
  int i;
  static int count_file=1;
  double x1[Npm],y1[Npm];
  char filename[128];   
  sprintf(filename,"%s/coord_%.d.dat",argv[1],count_file);
  ofstream file;
  file.open(filename);
  for(i=0;i<Np;i++){
    x1[i]=x[i]-x_corr;
    y1[i]=y[i]-y_corr;
  }
  p_bound(x1, y1, Np, L);  
  for(i=0;i<Np;i++)
    file << x1[i] << " " << y1[i]<< " "<<theta[i]<<" " << r1[i] <<" " << phi6[i] <<" " << phi4[i] <<" "<<phi2[i]<<" "<< s[i] << endl; 
  file.close(); 

  sprintf(filename,"%s/gr_%.d.dat",argv[1],count_file);
  file.open(filename);
  for(i=1;i<1500;i++)
    file << 0.02*i << " " << gr[i]/(count-count_th)<< " " << g6[i]/(count-count_th)<< " " << g4[i]/(count-count_th)<< " " << gs[i]/(count-count_th)<< " " << cq[i]/(count-count_th)<< endl; 
  file.close(); 
  count_file++;
  return 0;
}

int output_t(double *x,double *y,double *theta,double avU, double avK,int Np,double t,double time_stamp,double x_corr,double y_corr,char *argv[])
{
  int i;
  char filename[64];
  char filename2[64];
  FILE *fp;
  FILE *fp2;
  
  sprintf(filename,"%s/time_coord.dat",argv[1]);
  fp=fopen(filename,"a+");
  
  for(i=0;i<Np;i++)
    fprintf(fp,"%f\t%f\t%f\t%f\n",t-time_stamp,x[i]-x_corr,y[i]-y_corr,theta[i]);
  fclose(fp);
  
  sprintf(filename2,"%s/time_energy.dat",argv[1]);
  fp2=fopen(filename2,"a+");
  fprintf(fp2,"%f\t%f\t%f\t%f\t%f\t%f\n",t-time_stamp, avK,avU,x_corr,y_corr,theta[0]);
  fclose(fp2);

  return 0;
}

double unif_rand(double left, double right)
{
  return left + (right - left) * rand() / RAND_MAX;
}
double gaussian_rand(void)
{
  static double iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if (iset == 0) {
    do {
      v1 = unif_rand(-1, 1);
      v2 = unif_rand(-1, 1);
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 0.50;
    return v2 * fac;
  }
  else {
    iset = 0;
    return gset;
  }
}
int ini_coord_rand(double* x, double* y, double* a, int Np, double L, double r1, double r2,double* theta){
  int k, p;
 p = 0;
 srandom(time(0));

  for (k = 0; k < Np/2; k++) {
    x[2*k] = unif_rand(0, 1) * L;
    y[2*k] = unif_rand(0, 1) * L;
    theta[2*k] = unif_rand(0, 1) * pi;
    a[2*k] = r1;
    x[2*k+1] = unif_rand(0, 1) * L;
    y[2*k+1] = unif_rand(0, 1) * L;
    theta[2*k+1] = unif_rand(0, 1) * pi;
    a[2*k+1] = r2;
  }
  return 0;
}
int ini(double* vx, double* vy, int Np, double* omega) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
    omega[j]=0.0;
  }
  return 0;
}


int f(int i,int M)
{
  int k;
  
  k=i;
  
  if(k<0)
    k+=M;
  if(k>=M)
    k-=M;
  
  return k;
}
void update(double L,int Np,double *x,double *y,int M,double RCHK,int (*list)[Pm])
{
  int i,j,k;
  int nx,ny;
  int l,m;
  double dx,dy,r;
  
  int (*map)[Npm]=new int[M*M][Npm];
  
  for(i=0;i<M;i++)
    for(j=0;j<M;j++)
     	map[i+M*j][0]=0;
  
  for(i=0;i<Np;i++){
    nx=f((int)(x[i]*M/L),M);
    ny=f((int)(y[i]*M/L),M);
    
    for(m=ny-1;m<=ny+1;m++){
      for(l=nx-1;l<=nx+1;l++){
	map[f(l,M)+M*f(m,M)][map[f(l,M)+M*f(m,M)][0] +1]=i;
	map[f(l,M)+M*f(m,M)][0]++;
	//	printf("%d\n", map[f(l,M)+M*f(m,M)][0]);
      }
    }
  }
  
  for(i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)(x[i]*M/L),M);
    ny = f((int)(y[i]*M/L),M);
    
    for (k=1; k<=(map[nx+M*ny][0]); k++){
      j = map[nx+M*ny][k];
      if(j>i){
	dx =x[i] - x[j];
	dy =y[i] - y[j];
		
	if(dx<-L/2.0)
	  dx+=L;
	else if(dx> L/2.0)
	  dx-=L;

	if(dy<-L/2.0){
	  dy+=L;
	}
	else if(dy> L/2.0)
	  dy-=L;

	r = dx*dx + dy*dy;
	
	if(r<RCHK*RCHK){
	  list[i][0]++;
	  list[i][list[i][0]]=j;
	}
      }
    }
  }
  delete []map;
}


int calc_force_hs(double* x, double* y, double L, int Np, double* a, double* kx, double* ky, double *kth,double* avU,int (*list)[Pm],double *theta,double chi) {
  int i, j, k;
  *avU = 0.0;
  double r,r2,r4;
  double t, drU,U,Ai,Aj,A,sin_ij,cos_ij;
  double dx, dy;
  double aij;
  double cut;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    kth[k] = 0.0;
  }
  for (i = 0; i < Np; i++){
    for (j = 1; j <=list[i][0]; j++){
      dx = x[list[i][j]] - x[i];
      dy = y[list[i][j]] - y[i];
      if (dx > (0.5 * L))
	dx -= L;
      if (dx < -(0.5 * L))
	dx += L; 
      if (dy > (0.5 * L))
	dy -= L;
      if (dy < -(0.5 * L))
	dy += L;
      aij = (a[i] + a[list[i][j]]) / 2.0;
      r2 = dx * dx + dy * dy;
      r4=r2*r2;
 
      Ai=cos(theta[i])*dx+sin(theta[i])*dy;
      Aj=cos(theta[list[i][j]])*dx+sin(theta[list[i][j]])*dy;
      A= 1+chi*(Ai*Ai+Aj*Aj)/r2;
      cut = A*aij;
      if (r2 < cut*cut) {
	r=sqrt(r2);
	t = r / aij;

	drU = -(1 - t/A) /aij/A; //analytical calculation of the 1'st derivative
	U = 0.5*(1-t/A)*(1-t/A);

	cos_ij=cos(theta[i])+cos(theta[list[i][j]]);
	sin_ij=sin(theta[i])+sin(theta[list[i][j]]);   
	 
      }
      else {
	drU = 0.0;
	continue;
      } 
       
      kx[list[i][j]] -= drU * dx / r;
      kx[i] += drU * dx / r;
      ky[list[i][j]] -= drU * dy / r;
      ky[i] += drU * dy / r;
       


      kx[i]          += -drU*r/A*(
				  + 2.*chi*Ai*(cos(theta[i])*dy*dy-sin(theta[i])*dx*dy)/r4
				  +2.*chi*Aj*(cos(theta[list[i][j]])*dy*dy-sin(theta[list[i][j]])*dx*dy)/r4);
       
      kx[list[i][j]] -= -drU*r/A*(
				  +2.*chi*Ai*(cos(theta[i])*dy*dy-sin(theta[i])*dx*dy)/r4
				  +2.*chi*Aj*(cos(theta[list[i][j]])*dy*dy-sin(theta[list[i][j]])*dx*dy)/r4);
       
      ky[i]          += -drU*r/A*(
				  +2.*chi*Ai*(-cos(theta[i])*dx*dy+sin(theta[i])*dx*dx)/r4
				  + 2.*chi*Aj*(-cos(theta[list[i][j]])*dx*dy+sin(theta[list[i][j]])*dx*dx)/r4);
       
      ky[list[i][j]]  -= -drU*r/A*( 
				   +2.*chi*Ai*(-cos(theta[i])*dx*dy+sin(theta[i])*dx*dx)/r4
				   + 2.*chi*Aj*(-cos(theta[list[i][j]])*dx*dy+sin(theta[list[i][j]])*dx*dx)/r4);
       
       
      kth[i]         -= -t/A*drU*chi*(sin(2.*theta[i])*(-dx*dx+dy*dy)+2.*cos(2.*theta[i])*dx*dy)/r2;
      kth[list[i][j]]-= -t/A*drU*chi*(sin(2.*theta[list[i][j]])*(-dx*dx+dy*dy)+2.*cos(2.*theta[list[i][j]])*dx*dy)/r2;
       
      *avU += U;      
    }
  }
  *avU /= double(Np);
  return 0;
}

int calc_force(double* x, double* y, double L, int Np, double* a, double* kx, double* ky,double *kth,double* avU,double *avU2,int (*list)[Pm],double *theta,double chi) {
  int i, j, k;
  // *avU = 0.0;
  double r2,r4;
  double w2,w6,w12,drU,U12,A,Ai,Aj,cos_ij,sin_ij,fx_ij,fy_ij,U1;
  double dx, dy;
  double aij;
  double cut;
  // cut = 2.0;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    kth[k] = 0.0;
  }
 for (i = 0; i < Np; i++)
    {
      for (j = 1; j <=list[i][0]; j++)
	{
          dx = x[list[i][j]] - x[i];
          dy = y[list[i][j]] - y[i];
          if (dx > (0.5 * L))
            dx -= L;
          if (dx < -(0.5 * L))
            dx += L;
          if (dy > (0.5 * L))
            dy -= L;
          if (dy < -(0.5 * L))
            dy += L;
          aij = (a[list[i][j]] + a[i]) / 2.0;
          r2 = dx * dx + dy * dy; //avoid using sqrt() 
	  r4=r2*r2; 
          w2 = aij*aij / r2;
	  w6=w2*w2*w2;
	  w12=w6*w6;

	  Ai=cos(theta[i])*dx+sin(theta[i])*dy;
	  Aj=cos(theta[list[i][j]])*dx+sin(theta[list[i][j]])*dy;
	  A=  chi*(Ai*Ai+Aj*Aj)/r2;
	  drU = 4.*((-12.0) *(1.+A)*w12 + 6.0*w6)/ r2; //analytical calculation of the 1'st derivative / r
	  U12 = 4.*w12;
	
	  fx_ij=drU * dx
	    +2.*chi*Ai*(cos(theta[i])*dy*dy-sin(theta[i])*dx*dy)/r4*U12
	    +2.*chi*Aj*(cos(theta[list[i][j]])*dy*dy-sin(theta[list[i][j]])*dx*dy)/r4*U12; 
	  	  
	  fy_ij=drU * dy
	    +2.*chi*Ai*(-cos(theta[i])*dx*dy+sin(theta[i])*dx*dx)/r4*U12
	    +2.*chi*Aj*(-cos(theta[list[i][j]])*dx*dy+sin(theta[list[i][j]])*dx*dx)/r4*U12;
	  
	  
          if(r2< 3.0*3.0){
	    kx[i]+=fx_ij;
	    kx[list[i][j]]-=fx_ij;
	    	    
	    ky[i]+=fy_ij;	    
	    ky[list[i][j]]-=fy_ij;
	    
	    kth[i]         -=chi*(sin(2.*theta[i])*(-dx*dx+dy*dy)+2.*cos(2.*theta[i])*dx*dy)*U12/r2;
	    kth[list[i][j]]-=chi*(sin(2.*theta[list[i][j]])*(-dx*dx+dy*dy)+2.*cos(2.*theta[list[i][j]])*dx*dy)*U12/r2;
	    U1=4.*((1.+A)*w12-w6);
	    *avU +=U1/Np;
	    *avU2+=U1*U1/Np;      
	  }
        }
    }

 return 0; 
}

int eq_motion(double* x, double* y,double* theta, double* vx, double* vy, double *omega, double dt, double* kx, double* ky, double*kth, int Np, double* avK, double Th) {
  double zeta;
  zeta = 1.;
  int k;
  *avK=0.0;
  for (k = 0; k < Np; k++) {
    vx[k] += -vx[k] * zeta * dt + kx[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vy[k] += -vy[k] * zeta * dt + ky[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    omega[k]+= -omega[k]*zeta*dt+ kth[k]*dt+sqrt(2. * zeta * Th * dt) * gaussian_rand();
  
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    theta[k] +=omega[k] * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k] + omega[k] *omega[k];
  }
  *avK = *avK / Np / 2.0/1.5;
  return 0;
}

int main(int argc, char *argv[])
{
  srand((unsigned) time(NULL));
  double t,avU=0.0,avU2=0.0,avU0=0.0,avK=0.0,avK0=0.0,avphi6=0.0,avphi4=0.0,avphi2=0.0,avs=0.0,avphiq=0.0,dummy;
  double x_corr=0.0,y_corr=0.0;
  int i,count=0,count_num_update=0;
  double sampling_time,time_stamp=0.;
  double disp_max=0.0;
  double disp_ave=0.0;
  int Np = 16000;int count_th=200;
  double r1=1.0, r2=1.0;
  double disp_th =10.0;
  double dt=0.01, time_max = 1.e+6; //parameters
  double sampling_time_max=30000.;
  double time_coord= 1000;
  double time_stable_1 = 100.;
  double time_stable_2 = 1000.;
  double Th;
  double chi= atof(argv[1]);
  double timer;

  double chi0=0.2;
  if(chi>6)
    dt=0.002;
 
  if(chi>30)
    dt=0.001;

  double RCHK=4.5;
  double L = sqrt(double(Np)*pi*(r1*r1+r2*r2)/8.*pow(2.,0.333333)*pow((1.+2.*chi),1./6.) / 0.95); 
  int    M=(int)(L/RCHK);
  cout << "L=" << L <<" "<< "M="<<M <<endl;
  
  double* x, * y,* x0, * y0, * vx, * vy, * a, * kx, * ky,*kth,*x_update,*y_update,*theta,*omega,*phi6,*phi4,*phi2,*s,*gr,*g6,*g4,*g2,*gs,*cq,*qx,*qy;
  int (*list)[Pm]=new int[Npm][Pm];
  
  x = new double[Np];
  y = new double[Np];
  x0 = new double[Np];
  y0 = new double[Np];
  x_update = new double[Np];
  y_update = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  a = new double[Np];
  kx = new double[Np];
  ky = new double[Np];
  kth = new double[Np];
  phi6 = new double[Np];
  phi4 = new double[Np];
  qx = new double[500000];
  qy = new double[500000];
  phi2 = new double[Np];
  s= new double[Np];
  
  theta = new double[Np];
  omega = new double[Np];

  gr = new double[NR];
  g6 = new double[NR];
  gs = new double[NR];
  g4 = new double[NR];
  g2 = new double[NR];
  cq = new double[NR];

  char filename[128];
  
  ini_coord_rand(x, y, a, Np, L, r1, r2,theta);
  ini(vx, vy, Np,omega);  
  
  sprintf(filename,"%s/energy_time.txt",argv[1]);
  ofstream file;
  file.open(filename);
  
  //HP//////////////////////////////////////
  //  double chi0=0.2;
  double RCHK1=1.5;
  int M1=int(L/RCHK1);  
  for (t = 0.; t < time_stable_1; t += dt) {
    update(L,Np,x,y,M1,RCHK1,list);
    calc_force_hs(x, y, L, Np, a, kx, ky,kth, &dummy,list,theta,chi0);
    eq_motion(x, y, theta, vx, vy, omega, dt, kx, ky,kth, Np, &avK, 0.0);
    p_bound(x, y, Np, L);
    //cout << t <<endl;
  }

  ///////////////////////////////////////////

  update(L,Np,x,y,M,RCHK,list);
  avU0=0.0;

  for(Th=2.5;Th>=0.05;Th-=0.05){
    count=0;
    timer=0.0;
    avU=0.0,avU2=0.0,avphi6=0.0;avphi4=0.0;avphi2=0.0;avs=0.0;disp_ave=0.0;avK=0.0;avphiq=0.0;
    
    for(i=0;i<3000;i++){
      gr[i]=0.0;
      g6[i]=0.0;
      g4[i]=0.0;
      g2[i]=0.0;
      gs[i]=0.0;
      cq[i]=0.0;
    }
    
    copy(x_update,y_update,x,y,Np,x_corr,y_corr);
    copy(x0,y0,x,y,Np,x_corr,y_corr);
    
    for(;;){
      calc_force(x, y, L, Np, a, kx, ky,kth, &dummy,&dummy,list,theta,chi);
      eq_motion(x, y, theta, vx, vy, omega, dt, kx, ky,kth, Np, &avK0,Th);
      com_correction(x,y,&x_corr,&y_corr,Np, L);
      p_bound(x, y, Np, L);
      
      if(disp_ave> 0.001*0.001*(count+1)*(count+1)*disp_th*disp_th){
	//stage gate
	if(count==count_th){
	  sq(x,y,qx,qy,L,Np);
	  cout <<"sq OK Th="<<Th<<endl;
	}
	//production
	if(count>=count_th){
	  calc_order(x,y,theta,Np,L,phi6,phi4,phi2,s,&avphi6,&avphi4,&avphi2,&avs,a,gr,g6,g4,g2,gs,chi,qx,qy,cq,&avphiq);
	  calc_force(x, y, L, Np, a, kx, ky,kth, &avU,&avU2,list,theta,chi);
	  // cout <<count <<endl;
	}
	count++;
	avK+=avK0;
      }
      
      if(disp_ave>disp_th*disp_th){
	output(x,y,x_corr,y_corr,L,theta,a,Np,t,dt,phi6,phi4,phi2,s,argv,gr,g6,g4,gs,cq,count,avs,count_th);
	file <<Th<<" "<<avK/count<<" "<<avU/(count-count_th)<<" "<<(avU2/(count-count_th)-avU*avU/(count-count_th)/(count-count_th))/Th/Th<<" "<<avphi6/(count-count_th)<<" "<<avphi4/(count-count_th)<<" "<<avphi2/(count-count_th)<<" "<<avs/(count-count_th)<<" "<<avphiq/(count-count_th)<<endl;
	cout <<"T="<< Th << endl;
	avU0=avU;
	if(avphi6/(count-count_th)>0.2||avs/(count-count_th)>0.5)
	  disp_th=0.8;
	break;
      }
      
      ///////auto update////////////////////
      calc_disp(&disp_max,&disp_ave,x,y,x_update,y_update,x0,y0,Np,L,x_corr,y_corr);
      count_num_update++;
      if(disp_max>0.5*0.5){
	update(L,Np,x,y,M,RCHK,list);
	disp_max=0.0;
	copy(x_update,y_update,x,y,Np,x_corr,y_corr);
	count_num_update=0;
      }
      ////////////////////////////
    }
  }
  file.close();
  
  delete[] x;
  delete[] y;
  delete[] x0;
  delete[] y0;
  delete[] x_update;
  delete[] y_update;
  delete[] vx;
  delete[] vy;
  delete[] a;
  delete[] kx;
  delete[] ky;
  delete[] phi6;
  delete[] phi4;
  delete[] phi2;
  delete[] qx;
  delete[] qy;
  delete[] s;
  
  delete[] g6;
  delete[] g4;
  delete[] g2;
  delete[] cq;
  
  delete[] kth;
  delete[] list;
  delete[] theta;
  delete[] omega;
  return 0;
}
