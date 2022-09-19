# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <time.h>
#include <iomanip> // setprecisionを使用するのに必要
using namespace std;
# define Npm 16000
# define NR 20000
# define Nn 50
# define Pm 500
# define pi 3.141592653589793238462643383279
int copy(double *x_update,double *y_update,double *z_update,double *x,double *y,double *z,int Np,double x_corr,double y_corr,double z_corr){
 int i;
  for(i=0;i<Np;i++){
    x_update[i]=x[i];
    y_update[i]=y[i];
    z_update[i]=z[i];
  }
  return 0;
}

int calc_disp(double *disp_max,double *disp_ave,double *x,double *y,double *z,double *x_update,double *y_update,double *z_update,double *x0,double *y0,double *z0,int Np,int L,double x_corr,double y_corr,double z_corr){
  int i;
  double dx,dy,dz;
  double disp;
  
  *disp_ave = 0.0;

  for(i=0;i<Np;i++){
    dx=x[i]-x_update[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    dy=y[i]-y_update[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    dz=z[i]-z_update[i];
  
    disp = dx*dx + dy*dy + dz*dz;
    if(disp > *disp_max)
      *disp_max =disp;    
  }
  
  for(i=0;i<Np;i++){
    dx=x[i]-x0[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    dy=y[i]-y0[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    dz=z[i]-z0[i];
    disp = dx*dx + dy*dy + dz*dz;
    *disp_ave+=disp/Np;
  }  
  return 0;
}

int calc_disp_max(double *disp_max,double *x,double *y,double *z,double *x_update,double *y_update,double *z_update,int Np,int L){
  int i;
  double dx,dy,dz;
  double disp;
  
  for(i=0;i<Np;i++){
    dx=x[i]-x_update[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    
    dy=y[i]-y_update[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    
    dz=z[i]-z_update[i];
    
    disp = dx*dx+dy*dy+dz*dz;    
    if(disp > *disp_max)
      *disp_max =disp;    
  }
  
  return 0;
}


int p_bound(double* x, double* y, double* z,int Np, double L) {
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




int output(double *x,double *y,double *z,double* nx, double* ny,double* nz,double x_corr,double y_corr,double z_corr,double L,double *theta1,double *theta2,double *r1,int Np,double temp, double c,double kappa,double phi,double anc){
  int i;
  static int count_file=1;
  double x1[Npm],y1[Npm],z1[Npm];
  char filename[128];   
  sprintf(filename,"coord_anc%.1f_mix_phi%.2f_c%.1f_T%.1f_k%.1f_N%d.dat",anc,phi,c,temp,kappa,Np);
  ofstream file;
  file.open(filename);
  for(i=0;i<Np;i++){
    x1[i]=x[i];
    y1[i]=y[i];
    z1[i]=z[i];
    x1[i]-=L*floor(x1[i]/L);
    y1[i]-=L*floor(y1[i]/L);
    z1[i]-=L*floor(z1[i]/L);
  }
  for(i=0;i<Np;i++)
    file <<r1[i]<<" "<< x1[i] << " " << y1[i] << " " << z1[i]<< " "<<nx[i]<<" " <<ny[i]<<" " << nz[i] << endl; 
  file.close(); 
  
  count_file++;
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
      rsq = v1 * v1 + v2 * v2 ;
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
void affine_coord(double* x,double* y,double* z,double phi,int Np,double dphi){
  for(int k = 0; k < Np; k++){
    x[k] = pow(phi/(phi + dphi),1./2.)*x[k];
    y[k] = pow(phi/(phi + dphi),1./2.)*y[k];
  }
}


void ini_hex(double *x,double *y,double*z, double* nx, double* ny,double*nz,double *theta,double *a,double L,int Np,double kappa){
  int num_x = (int)sqrt(Np)+1;
  int num_y = (int)sqrt(Np)+1;
  
  double lx = L/double(num_x);
  double ly = L/double(num_y);
  double n;
  int i,j,k=0;
  double shift;
  for(j=0;j<num_y;j++){
    for(i=0;i<num_x;i++){
      x[i+num_x*j] = lx*i + lx/2.;
      y[i+num_x*j] = ly*j + ly/2.;
      z[i+num_x*j] = kappa*pow(2.,1./6.)/2.;
      nx[i+num_x*j] = unif_rand(-1., 1.);
      ny[i+num_x*j] = unif_rand(-1., 1.);
      nz[i+num_x*j] = unif_rand(-1., 1.);
      a[i+num_x*j]=1.0;
      k++;
      if(k>=Np)
        break;
    }
    if(k>=Np)
      break;
  }
}  

int initial(double* nx, double* ny,double*nz,int Np){
  for(int k = 0; k <Np;k++){
    double n = sqrt(nx[k] * nx[k] + ny[k] * ny[k] + nz[k] * nz[k]);
    nx[k] = nx[k] / n;
    ny[k] = ny[k] / n;
    nz[k] = nz[k] / n;
  }
  return 0;
}

int ini(double* vx, double* vy,double* vz, double* nx_pr, double* ny_pr,double* nz_pr, int Np, double* omega) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
    vz[j] = 0.0;
    nx_pr[j] = 0.0;
    ny_pr[j] = 0.0;
    nz_pr[j] = 0.0;
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


void update(double L,double Lz,int Np,double *x,double *y,double *z,int M,int Mz,double Rcell,int (* list)[Pm]) 
{
  int i,j,k;
  int nx,ny,nz;
  int l,m,n;
  double dx,dy,dz,r;
  
  int (*map)[Npm]=new int[M*M*Mz][Npm];
  for(k=0;k<Mz;k++)
    for(j=0;j<M;j++)
      for(i=0;i<M;i++)
	map[i+M*j+M*M*k][0]=0;
  
  for(i=0;i<Np;i++){
    nx=f((int)(x[i]*M/L),M);
    ny=f((int)(y[i]*M/L),M);
    nz=f((int)(z[i]*Mz/L),Mz); //address of particle i is determined. 
    for(n=nz-1;n<=nz+1;n++)
      for(m=ny-1;m<=ny+1;m++)
	for(l=nx-1;l<=nx+1;l++){
	  map[f(l,M)+M*f(m,M)+M*M*f(n,Mz)][map[f(l,M)+M*f(m,M)+M*M*f(n,Mz)][0] +1]=i; 
	  map[f(l,M)+M*f(m,M)+M*M*f(n,Mz)][0]++;
        }
  } 
  
    
  for (i=0;i<Np;i++){ 
    list[i][0]=0;
    nx = f((int)(x[i]*M/L),M); 
    ny = f((int)(y[i]*M/L),M); 
    nz = f((int)(z[i]*Mz/L),Mz);
    
    for (k=1; k<=(map[nx+M*ny+M*M*nz][0]); k++){ 
      j = map[nx+M*ny+M*M*nz][k];
      if(j>i) {
	dx =x[i] - x[j]; 
	dy =y[i] - y[j]; 
	dz =z[i] - z[j];
	
	if(dx<-L/2.0) 
	  dx+=L;
	else if(dx> L/2.0) 
	  dx-=L;
	if(dy<-L/2.0)
	  dy+=L;
	else if(dy> L/2.0)
	  dy-=L;
	
	r = dx*dx + dy*dy + dz*dz;
	
	if(r<Rcell*Rcell){ 
	  list[i][0]++;
	  list[i][list[i][0]]=j;
	}
      }
    }
  }
  delete
    []map;
}

int calc_force(double* x, double* y, double* z, double L, double* nx, double* ny, double* nz,int Np, double* a, double* kx, double* ky,double* kz,double *knx,double *kny,double *knz,double*kw,double* kx_w, double* ky_w,double* kz_w,double *knx_w,double *kny_w,double *knz_w,int (*list)[Pm],double *theta1,double *theta2,double chi,double chi_pr,double mu,double eta,double kappa,double* kx_c, double* ky_c,double* kz_c,double *knx_c,double *kny_c,double *knz_c,double c,double*p,double* avU,double *avU2,double anc) {
  int i, j, k;
  // *avU = 0.0;
  double r,r2,r3,r4,ri,rj,cij,dij;
  double R2p,R2n,R1p,R1n,R2p_pr,R2n_pr,R1p_pr,R1n_pr;
  double t,drU,diU,djU,dcU,U,A,B,W,X,A1,A2,A6,A12,W2,W6,W12,W13,X2,X6,X12,X13,fx_ij,fy_ij,fz_ij,fw_i,fx_ij_c,fy_ij_c,fz_ij_c,U1,Lz,A_c,B_c,drUc,diUc,djUc,dcUc,ddUc,Uc;
  double e1,e2;
  double dx, dy, dz;
  double aij, aij_pr, aij3, aij_pr3;
  double cut;

  double rw_low,r2w_low,r3w_low,r4w_low,riw_low,rjw_low,cijw,rjw_high,cijw_high;
  double rw_high,r2w_high,r3w_high,r4w_high,riw_high;
  double R2pw_low,R2nw_low,R1pw_low,R1nw_low, R2p_pr_w_low, R2n_pr_w_low, R1p_pr_w_low, R1n_pr_w_low;
  double R2pw_high,R2nw_high,R1pw_high,R1nw_high, R2p_pr_w_high, R2n_pr_w_high, R1p_pr_w_high, R1n_pr_w_high;
  double drU_w_low,diU_w_low,djU_w_low,drU_w_high,diU_w_high,djU_w_high,dcU_w_low,dcU_w_high,Aw_low,Bw_low,A1w_low,A2w_low,A6w_low,A12w_low,Aw_high,Bw_high,A1w_high,A2w_high,A6w_high,A12w_high,fx_ij_w_low,fy_ij_w_low,fz_ij_w_low,fx_ij_w_high,fy_ij_w_high,fz_ij_w_high;
  double e1w,e2w_low,e2w_high;
  double  lw_low, lw_high; 
  double aij_w, aij_pr_w_low, aij_pr_w_high, aij_w3, aij_pr3_w_low, aij_pr3_w_high;

  Lz = kappa*pow(2.,1./6.);
  double V = L*L*Lz;

  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    kz[k] = 0.0;
    knx[k] = 0.0;
    kny[k] = 0.0;
    knz[k] = 0.0;
    kw[k] = 0.0;
    kx_c[k] = 0.0;
    ky_c[k] = 0.0;
    kz_c[k] = 0.0;
    knx_c[k] = 0.0;
    kny_c[k] = 0.0;
    knz_c[k] = 0.0;
  }

    for (k = 0; k < Np; k++) {
      kx_w[k] = 0.0;
      ky_w[k] = 0.0;
      kz_w[k] = 0.0;
      knx_w[k] = 0.0;
      kny_w[k] = 0.0;
      knz_w[k] = 0.0;
    }
    for (i = 0; i < Np; i++)
      {
	lw_low = z[i] + kappa/2.*pow(2.,1./6.);
	lw_high = z[i] - Lz - kappa/2.*pow(2.,1./6.);
	aij_w = 1.0;
	aij_w3 = aij_w*aij_w*aij_w;
	
	// cout<<theta[list[i][j]]<<endl;
	rw_low = sqrt(lw_low*lw_low);
	r2w_low = lw_low*lw_low;
	r3w_low = r2w_low*rw_low;
	r4w_low = r2w_low*r2w_low;
	
	rw_high = sqrt(lw_high*lw_high);
	r2w_high = lw_high*lw_high;
	r3w_high = r2w_high*rw_high;
	r4w_high = r2w_high*r2w_high;
	
	// n of j is always (0,0,1)//  MIRR0R ANCHORING
	
	riw_low =  lw_low*nz[i];
	rjw_low = lw_low;
	cijw = nz[i];
	
	riw_high =  lw_high*nz[i];
	rjw_high = lw_high;
	
	R2pw_low =(riw_low + rjw_low)*(riw_low + rjw_low)/(1 + chi*cijw);
	R2nw_low =(riw_low - rjw_low)*(riw_low - rjw_low)/(1 - chi*cijw);
	R1pw_low =(riw_low + rjw_low)/(1 + chi*cijw);
	R1nw_low =(riw_low - rjw_low)/(1 - chi*cijw);
	
	R2p_pr_w_low =(riw_low + rjw_low)*(riw_low + rjw_low)/(1 + chi_pr*cijw);
	R2n_pr_w_low =(riw_low - rjw_low)*(riw_low - rjw_low)/(1 - chi_pr*cijw);
	R1p_pr_w_low =(riw_low + rjw_low)/(1 + chi_pr*cijw);
	R1n_pr_w_low =(riw_low - rjw_low)/(1 - chi_pr*cijw); 
	
	aij_pr_w_low = aij_w/sqrt(1.-(chi/2./r2w_low)*(R2pw_low+R2nw_low));
	aij_pr3_w_low = aij_pr_w_low* aij_pr_w_low* aij_pr_w_low;
	
	R2pw_high =(riw_high + rjw_high)*(riw_high + rjw_high)/(1 + chi*cijw);
	R2nw_high =(riw_high - rjw_high)*(riw_high - rjw_high)/(1 - chi*cijw);
	R1pw_high =(riw_high + rjw_high)/(1 + chi*cijw);
	R1nw_high =(riw_high - rjw_high)/(1 - chi*cijw);
	
	R2p_pr_w_high =(riw_high + rjw_high)*(riw_high + rjw_high)/(1 + chi_pr*cijw);
	R2n_pr_w_high =(riw_high - rjw_high)*(riw_high - rjw_high)/(1 - chi_pr*cijw);
	R1p_pr_w_high =(riw_high + rjw_high)/(1 + chi_pr*cijw);
	R1n_pr_w_high =(riw_high - rjw_high)/(1 - chi_pr*cijw);
	
	aij_pr_w_high = aij_w/sqrt(1.-(chi/2./r2w_high)*(R2pw_high+R2nw_high));
	aij_pr3_w_high = aij_pr_w_high*aij_pr_w_high*aij_pr_w_high;
	
	A1w_low = (aij_w/(rw_low-aij_pr_w_low+aij_w));
	A2w_low =A1w_low*A1w_low;
	A6w_low = A2w_low*A2w_low*A2w_low;
	A12w_low = A6w_low*A6w_low;
	
	A1w_high = (aij_w/(rw_high-aij_pr_w_high+aij_w));
	A2w_high =A1w_high*A1w_high;
	A6w_high = A2w_high*A2w_high*A2w_high;
	A12w_high = A6w_high*A6w_high;
	
	Aw_low = A12w_low - A6w_low;
	Bw_low = 12.*A12w_low*A1w_low - 6.*A6w_low*A1w_low;
	
	Aw_high = A12w_high - A6w_high;
	Bw_high = 12.*A12w_high*A1w_high - 6.*A6w_high*A1w_high;
	
	//	    e1 = 1/sqrt(1-chi*chi*cij*cij);  //epsilon1
	e1w = 1/sqrt(1-chi*chi*cijw*cijw);  //epsilon1
	e2w_low = 1.-(chi_pr/2./r2w_low)*(R2p_pr_w_low + R2n_pr_w_low); // epsilon2
	e2w_high = 1.-(chi_pr/2./r2w_high)*(R2p_pr_w_high + R2n_pr_w_high); // epsilon2
	
	
	drU_w_low = 4.*pow(e1w,eta)*pow(e2w_low,mu)*(((Aw_low*mu*chi_pr)/e2w_low/r3w_low)*(R2p_pr_w_low + R2n_pr_w_low)-((aij_pr3_w_low*chi*Bw_low)/2./aij_w3/r3w_low)*(R2pw_low + R2nw_low)-Bw_low/aij_w) ;
	diU_w_low = 4.*pow(e1w,eta)*pow(e2w_low,mu)*(-((Aw_low*mu*chi_pr)/e2w_low/r2w_low)*(R1p_pr_w_low + R1n_pr_w_low)+((aij_pr3_w_low*chi*Bw_low)/2./aij_w3/r2w_low)*(R1pw_low+R1nw_low));
	djU_w_low = 4.*pow(e1w,eta)*pow(e2w_low,mu)*(-((Aw_low*mu*chi_pr)/e2w_low/r2w_low)*(R1p_pr_w_low - R1n_pr_w_low)+((aij_pr3_w_low*chi*Bw_low)/2./aij_w3/r2w_low)*(R1pw_low-R1nw_low));//analytical calculation of the 1'st derivative / r*nj
	
	drU_w_high = 4.*pow(e1w,eta)*pow(e2w_high,mu)*(((Aw_high*mu*chi_pr)/e2w_high/r3w_high)*(R2p_pr_w_high + R2n_pr_w_high)-((aij_pr3_w_high*chi*Bw_high)/2./aij_w3/r3w_high)*(R2pw_high + R2nw_high)-Bw_high/aij_w) ;//analytical calculation of the 1'st derivative / r
	diU_w_high = 4.*pow(e1w,eta)*pow(e2w_high,mu)*(-((Aw_high*mu*chi_pr)/e2w_high/r2w_high)*(R1p_pr_w_high + R1n_pr_w_high)+((aij_pr3_w_high*chi*Bw_high)/2./aij_w3/r2w_high)*(R1pw_high+R1nw_high));//analytical calculation of the 1'st derivative / r*ni
	djU_w_high = 4.*pow(e1w,eta)*pow(e2w_high,mu)*(-((Aw_high*mu*chi_pr)/e2w_high/r2w_high)*(R1p_pr_w_high - R1n_pr_w_high)+((aij_pr3_w_high*chi*Bw_high)/2./aij_w3/r2w_high)*(R1pw_high-R1nw_high));//analytical calculation of the 1'st derivative / r*nj
	
	dcU_w_low= 4.*pow(e1w,eta)*pow(e2w_low,mu)*(((Aw_low*mu*chi_pr*chi_pr)/2./e2w_low/r2w_low)*(R1p_pr_w_low*R1p_pr_w_low-R1n_pr_w_low*R1n_pr_w_low)-((aij_pr3_w_low*chi*chi*Bw_low)/4./aij_w3/r2w_low)*(R1pw_low*R1pw_low-R1nw_low*R1nw_low)+Aw_low*eta*chi*chi*cijw*e1w*e1w);//analytical calculation of the 1'st derivative / ri*ni
	dcU_w_high= 4.*pow(e1w,eta)*pow(e2w_high,mu)*(((Aw_high*mu*chi_pr*chi_pr)/2./e2w_high/r2w_high)*(R1p_pr_w_high*R1p_pr_w_high-R1n_pr_w_high*R1n_pr_w_high)-((aij_pr3_w_high*chi*chi*Bw_high)/4./aij_w3/r2w_high)*(R1pw_high*R1pw_high-R1nw_high*R1nw_high)+Aw_high*eta*chi*chi*cijw*e1w*e1w);//analytical calculation of the 1'st derivative / ri*ni
	
	
	fx_ij_w_low = diU_w_low*nx[i];//+djU*nx[list[i][j]];  //du/dr
	fy_ij_w_low = diU_w_low*ny[i];//+djU*ny[list[i][j]];
	fz_ij_w_low= drU_w_low*lw_low/rw_low + diU_w_low*nz[i] + djU_w_low;//+djU*nz[list[i][j]];
	
	fx_ij_w_high = diU_w_high*nx[i];//+djU*nx[list[i][j]];  //du/dr
	fy_ij_w_high = diU_w_high*ny[i];//+djU*ny[list[i][j]];
	fz_ij_w_high = drU_w_high*lw_high/rw_high + diU_w_high*nz[i] + djU_w_high;//+djU*nz[list[i][j]];
	
	kx_w[i] = -anc*(fx_ij_w_low + fx_ij_w_high);
	ky_w[i] = -anc*(fy_ij_w_low + fy_ij_w_high);	    
	kz_w[i] = -anc*(fz_ij_w_low + fz_ij_w_high);
	
	knx_w[i] = anc*(diU_w_low*riw_low*nx[i] + dcU_w_low*cijw*nx[i] + diU_w_high*riw_high*nx[i] + dcU_w_high*cijw*nx[i]);
	kny_w[i] = anc*(diU_w_low*riw_low*ny[i] + dcU_w_low*cijw*ny[i] + diU_w_high*riw_high*ny[i] + dcU_w_high*cijw*ny[i]);
	knz_w[i] = anc*(-diU_w_low*(lw_low - riw_low*nz[i]) - dcU_w_low*(1.0 - cijw*nz[i]) - diU_w_high*(lw_high - riw_high*nz[i]) - dcU_w_high*(1.0 - cijw*nz[i]));
	
	for (j = 1; j <=list[i][0]; j++)
	  {  
	    dx = x[i] - x[list[i][j]];
	    dy = y[i] - y[list[i][j]];
	    dz = z[i] - z[list[i][j]];
	    if (dx > (0.5 * L))
	      dx -= L;
	    if (dx < -(0.5 * L))
	      dx += L;
	    if (dy > (0.5 * L))
	      dy -= L;
	    if (dy < -(0.5 * L))
	      dy += L;
	    if (dz > (0.5 * L))
	      dz -= L;
	    if (dz < -(0.5 * L))
	      dz += L;
	    
	    aij = (a[list[i][j]] + a[i]) / 2.0;
	    aij3 = aij*aij*aij;
	    
	    // cout<<theta[list[i][j]]<<endl;
	    r = sqrt(dx*dx + dy*dy + dz*dz);
	    r2 = dx*dx + dy*dy + dz*dz;
	    r3 = r2*r;
	    r4 = r2*r2;
	    
	    ri = dx*nx[i] + dy*ny[i] + dz*nz[i];
	    rj = dx*nx[list[i][j]] + dy*ny[list[i][j]] + dz*nz[list[i][j]];
	    cij = nx[i]*nx[list[i][j]] + ny[i]*ny[list[i][j]] + nz[i]*nz[list[i][j]];
	    
	    R2p =(ri+rj)*(ri+rj)/(1+chi*cij);
	    R2n =(ri-rj)*(ri-rj)/(1-chi*cij);
	    R1p =(ri+rj)/(1+chi*cij);
	    R1n =(ri-rj)/(1-chi*cij);
	    
	    aij_pr = aij/sqrt(1.-(chi/2./r2)*(R2p+R2n));
	    cut = pow(2.,1./6.)-1+aij_pr/aij;
	    
	    if(r2< cut*cut){    
	      dij = dx*(ny[i]*nz[list[i][j]] - nz[i]*ny[list[i][j]]) + dy*(nz[i]*nx[list[i][j]] - nx[i]*nz[list[i][j]]) + dz*(nx[i]*ny[list[i][j]] - ny[i]*nx[list[i][j]]);
	      
	      R2p_pr =(ri+rj)*(ri+rj)/(1+chi_pr*cij);
	      R2n_pr =(ri-rj)*(ri-rj)/(1-chi_pr*cij);
	      R1p_pr =(ri+rj)/(1+chi_pr*cij);
	      R1n_pr =(ri-rj)/(1-chi_pr*cij);    
	      
	      aij_pr3 = aij_pr*aij_pr*aij_pr;
	      
	      A1 = (aij/(r-aij_pr+aij));
	      A2 =A1*A1;
	      A6 = A2*A2*A2;
	      A12 = A6*A6;
	      
	      A = A12-A6;
	      B = 12.*A12*A1-6.*A6*A1;
	      
	      A_c = A1*A6;
	      B_c = 7.*A6*A2;
	      
	      //	    e1 = 1/sqrt(1-chi*chi-cij*cij);  //epsilon1
	      e1 = 1./sqrt(1-chi*chi*cij*cij);  //epsilon1
	      e2 = 1.-(chi_pr/2./r2)*(R2p_pr+R2n_pr); // epsilon2
	      
	      //start derivation
	      
	      drU = 4.*pow(e1,eta)*pow(e2,mu)*(((A*mu*chi_pr)/e2/r3)*(R2p_pr+R2n_pr)-((aij_pr3*chi*B)/2./aij3/r3)*(R2p+R2n)-B/aij);//analytical calculation of the 1'st derivative / r
	      diU = 4.*pow(e1,eta)*pow(e2,mu)*(-((A*mu*chi_pr)/e2/r2)*(R1p_pr+R1n_pr)+((aij_pr3*chi*B)/2./aij3/r2)*(R1p+R1n));//analytical calculation of the 1'st derivative / r*ni
	      djU = 4.*pow(e1,eta)*pow(e2,mu)*(-((A*mu*chi_pr)/e2/r2)*(R1p_pr-R1n_pr)+((aij_pr3*chi*B)/2./aij3/r2)*(R1p-R1n));//analytical calculation of the 1'st derivative / r*nj
	      dcU = 4.*pow(e1,eta)*pow(e2,mu)*(((A*mu*chi_pr*chi_pr)/2./e2/r2)*(R1p_pr*R1p_pr-R1n_pr*R1n_pr)-((aij_pr3*chi*chi*B)/4./aij3/r2)*(R1p*R1p-R1n*R1n)+A*eta*chi*chi*cij*e1*e1);//analytical calculation of the 1'st derivative / ri*ni
	      
	      fx_ij = drU*dx/r+diU*nx[i]+djU*nx[list[i][j]];  //du/dr
	      fy_ij = drU*dy/r+diU*ny[i]+djU*ny[list[i][j]];
	      fz_ij = drU*dz/r+diU*nz[i]+djU*nz[list[i][j]];
	      
	      kx[i] -= fx_ij;
	      kx[list[i][j]] += fx_ij;
	      
	      ky[i] -= fy_ij;	    
	      ky[list[i][j]] += fy_ij;
	      
	      kz[i] -= fz_ij;
	      kz[list[i][j]] += fz_ij;
	      
	    
	    
	      knx[i] -= diU*(dx - ri*nx[i]) + dcU*(nx[list[i][j]] - cij*nx[i]);
	      knx[list[i][j]] -= djU*(dx - rj*nx[list[i][j]]) + dcU*(nx[i] - cij*nx[list[i][j]]);
	      
	      kny[i] -= diU*(dy - ri*ny[i]) + dcU*(ny[list[i][j]] - cij*ny[i]);
	      kny[list[i][j]] -= djU*(dy - rj*ny[list[i][j]]) + dcU*(ny[i] - cij*ny[list[i][j]]);
	      
	      knz[i] -= diU*(dz - ri*nz[i]) + dcU*(nz[list[i][j]] - cij*nz[i]);
	      knz[list[i][j]] -= djU*(dz - rj*nz[list[i][j]]) + dcU*(nz[i] - cij*nz[list[i][j]]);
	      
	      // Chiral//
	      
	      drUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*cij*(((A_c*mu*chi_pr)/e2/r4)*(R2p_pr + R2n_pr) - ((aij_pr3*chi*B_c)/2./aij3/r4)*(R2p + R2n) - B_c/aij/r - A_c/r2);//analytical calculation of the 1'st derivative / r
	      diUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*cij*(-((A_c*mu*chi_pr)/e2/r3)*(R1p_pr + R1n_pr) + ((aij_pr3*chi*B_c)/2./aij3/r3)*(R1p + R1n));//analytical calculation of the 1'st derivative / r*ni
	      djUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*cij*(-((A_c*mu*chi_pr)/e2/r3)*(R1p_pr - R1n_pr) + ((aij_pr3*chi*B_c)/2./aij3/r3)*(R1p - R1n));//analytical calculation of the 1'st derivative / r*ni
	      dcUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*(((A_c*mu*chi_pr*chi_pr*cij)/2./e2/r3)*(R1p_pr*R1p_pr - R1n_pr*R1n_pr) - ((aij_pr3*chi*chi*B_c*cij)/4./aij3/r3)*(R1p*R1p - R1n*R1n) + (A_c*eta*chi*chi*cij*cij*e1*e1)/r);//analytical calculation of the 1'st derivative / ri*ni
	      ddUc = 4.*pow(e1,eta)*pow(e2,mu)*cij*(A_c/r);
	      
	      fx_ij_c = drUc*dx/r + diUc*nx[i] + djUc*nx[list[i][j]] + ddUc*(ny[i]*nz[list[i][j]] - nz[i]*ny[list[i][j]]);  //du/dr
	      fy_ij_c = drUc*dy/r + diUc*ny[i] + djUc*ny[list[i][j]] + ddUc*(nz[i]*nx[list[i][j]] - nx[i]*nz[list[i][j]]);
	      fz_ij_c = drUc*dz/r + diUc*nz[i] + djUc*nz[list[i][j]] + ddUc*(nx[i]*ny[list[i][j]] - ny[i]*nx[list[i][j]]);
	      
	      kx_c[i] -= fx_ij_c;
	      kx_c[list[i][j]] += fx_ij_c;
	      
	      ky_c[i] -= fy_ij_c;	    
	      ky_c[list[i][j]] += fy_ij_c;
	      
	      kz_c[i] -= fz_ij_c;
	      kz_c[list[i][j]] += fz_ij_c;
	      
	      
	      knx_c[i] -= diUc*(dx - ri*nx[i]) + dcUc*(nx[list[i][j]] - cij*nx[i]) + ddUc*(ri*(ny[list[i][j]]*nz[i] - nz[list[i][j]]*ny[i]) - cij*(dy*nz[i] - dz*ny[i]));
	      knx_c[list[i][j]] -= djUc*(dx - rj*nx[list[i][j]]) + dcUc*(nx[i] - cij*nx[list[i][j]]) + ddUc*(-rj*(-ny[list[i][j]]*nz[i] + nz[list[i][j]]*ny[i]) + cij*(dy*nz[list[i][j]] - dz*ny[list[i][j]]));
	      
	      kny_c[i] -= diUc*(dy - ri*ny[i]) + dcUc*(ny[list[i][j]] - cij*ny[i]) + ddUc*(ri*(nz[list[i][j]]*nx[i] - nx[list[i][j]]*nz[i]) - cij*(dz*nx[i] - dx*nz[i]));
	      kny_c[list[i][j]] -= djUc*(dy - rj*ny[list[i][j]]) + dcUc*(ny[i] - cij*ny[list[i][j]]) + ddUc*(-rj*(-nz[list[i][j]]*nx[i] + nx[list[i][j]]*nz[i]) + cij*(dz*nx[list[i][j]] - dx*nz[list[i][j]]));
	      
	      knz_c[i] -= diUc*(dz - ri*nz[i]) + dcUc*(nz[list[i][j]] - cij*nz[i]) + ddUc*(ri*(nx[list[i][j]]*ny[i] - ny[list[i][j]]*nx[i]) - cij*(dx*ny[i] - dy*nx[i]));
	      knz_c[list[i][j]] -= djUc*(dz - rj*nz[list[i][j]]) + dcUc*(nz[i] - cij*nz[list[i][j]]) + ddUc*(-rj*(-nx[list[i][j]]*ny[i] + ny[list[i][j]]*nx[i]) + cij*(dx*ny[list[i][j]] - dy*nx[list[i][j]]));
	      
	      *p += -(1./3.)*(dx*fx_ij + dy*fy_ij + dz*fz_ij + c*(dx*fx_ij_c + dy*fy_ij_c + dz*fz_ij_c));
	      
	      U1=4.*pow(e1,eta)*pow(e2,mu)*A + 4.*c*pow(e1,eta)*pow(e2,mu)*cij*dij*A_c/r;
	      *avU +=U1/Np;
	      *avU2+=U1*U1/Np;   
	    }
	  }
      }
    *p = *p/V;
    *avU /=Np;
	*avU2 /=Np;   
    
    return 0; 
}

int eq_motion(double* x, double* y,double* z,double* theta1,double* theta2, double* vx, double* vy, double* vz,double* nx,double* ny,double* nz, double* nx_pr,double* ny_pr,double* nz_pr, double* omega, double dt, double* kx, double* ky, double* kz,double* knx,double* kny,double* knz,double* kw ,double* kx_w, double* ky_w,double* kz_w,double *knx_w,double *kny_w,double *knz_w,int Np, double* avK, double Th,double* kx_c, double* ky_c, double* kz_c,double* knx_c,double* kny_c,double* knz_c, double c) {
  double zeta, zeta_r;
  double n;
  double n_pr2;
  double Rx;
  double Ry;
  double Rz;
  int k;
  
  zeta = 100.;
  zeta_r = 100.;
  *avK=0.0;
  for (k = 0; k < Np; k++) {
    vx[k] += -vx[k] * zeta * dt + kx[k] * dt + c * kx_c[k] * dt + kx_w[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vy[k] += -vy[k] * zeta * dt + ky[k] * dt + c * ky_c[k] * dt + ky_w[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vz[k] += -vz[k] * zeta * dt + kz[k] * dt + c * kz_c[k] * dt + kz_w[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    //    Rx[k] = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    //    Ry[k] = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    //    Rz[k] = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    Rx = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    Ry = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    Rz = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    n_pr2 = nx_pr[k]*nx_pr[k] + ny_pr[k]*ny_pr[k] + nz_pr[k]*nz_pr[k];
    // cout << "n_pr2=" << " " << n_pr2  <<  endl;
    //   nx_pr[k] += -n_pr2[k]*nx[k]*dt - zeta_r*nx_pr[k]*dt + knx[k]*dt + Ry[k]*nz[k] - Rz[k]*ny[k] ;
    //    ny_pr[k] += -n_pr2[k]*ny[k]*dt - zeta_r*ny_pr[k]*dt + kny[k]*dt + Rz[k]*nx[k] - Rx[k]*nz[k] ;
    //    nz_pr[k] += -n_pr2[k]*nz[k]*dt - zeta_r*nz_pr[k]*dt + knz[k]*dt + Rx[k]*ny[k] - Ry[k]*nx[k] ;
    nx_pr[k] += -n_pr2*nx[k]*dt - zeta_r*nx_pr[k]*dt + knx[k]*dt + c*knx_c[k]*dt + knx_w[k]*dt + Ry*nz[k] - Rz*ny[k] ;
    ny_pr[k] += -n_pr2*ny[k]*dt - zeta_r*ny_pr[k]*dt + kny[k]*dt + c*kny_c[k]*dt + kny_w[k]*dt + Rz*nx[k] - Rx*nz[k] ;
    nz_pr[k] += -n_pr2*nz[k]*dt - zeta_r*nz_pr[k]*dt + knz[k]*dt + c*knz_c[k]*dt + knz_w[k]*dt + Rx*ny[k] - Ry*nx[k] ;
    
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    z[k] += vz[k] * dt;
    
    nx[k] += nx_pr[k] * dt;
    ny[k] += ny_pr[k] * dt;
    nz[k] += nz_pr[k] * dt;
    
    n = sqrt(nx[k] * nx[k] + ny[k] * ny[k] + nz[k] * nz[k]);
    nx[k] = nx[k] / n;
    ny[k] = ny[k] / n;
    nz[k] = nz[k] / n;
    *avK += vx[k]*vx[k] + vy[k]*vy[k] + vz[k]*vz[k] + (nx_pr[k]*nx_pr[k] + ny_pr[k]*ny_pr[k] + nz_pr[k]*nz_pr[k])/(n*n);
  }
  *avK = *avK / Np / 2.0/ 2.5;
  // cout << *avK << endl;
  return 0;
}

int main(int argc, char *argv[])
{
  srand((unsigned) time(NULL));
  double t,avU=0.0,avU2=0.0,avU0=0.0,avK=0.0,avK0=0.0,p0=0.0,p=0.0,dummy;
  double x_corr=0.0,y_corr=0.0,z_corr=0.0;
  int i,count=0,count_num_update=0;
  
  double disp_max=0.0,disp_th2 =10.;
  double disp_ave;
  int Np = 2000;//10000;
  int count_th=200;
  //double disp_th2 = 180;
  double dt =0.0001;//  //parameters;
  double temp = 0.5;
  double Th;
  double phi=0.005;
  double dphi = 0.0005;
  double phi_final;
  
  double kappa = 4.4; // shape anisotropy parameter
  double kappa_pr = 0.5; // energy anisotropy parameter
  double chi = (kappa*kappa-1.)/(kappa*kappa+1.); //atof(argv[1]);  
  double c = 1.0;
  double anc = 0.0;
  
  double eta = 2.0;
  double mu = 1.0;
  double chi_pr = (pow(kappa_pr,1./mu)-1.)/(pow(kappa_pr,1./mu)+1.) ; 
  double timer;
  double RCHK=kappa+2.0;// cut off +1.0
  double L,Lz;
  L = sqrt(pow(2.,1./2.)*pi*kappa*Np/(6.*kappa*pow(2.,1./6.)*phi));
  //L = sqrt(pow(2.,1./2.)*pi*kappa*Np/(6.*kappa*phi));
  Lz = 20.;
  int M,Mz;
  M=(int)(L/RCHK);
  Mz=(int)(Lz/RCHK);
  
  double lz=kappa*pow(2.,1./6.);
  double rho=Np/L/L/lz;
  cout << "L=" << L <<" "<< "M="<<M <<endl;
  
  double* x, * y, *z,* x0, * y0, *z0,* vx, * vy,*vz, * a, * kx, * ky,*kz, * knx,* kny,* knz, * kx_w, * ky_w,*kz_w, * knx_w,* kny_w,* knz_w, * nx,* ny,* nz, * nx_pr,* ny_pr,* nz_pr,*x_update,*y_update,*z_update, *theta1, *theta2 ,*omega,*r1,*kw,* kx_c, * ky_c,*kz_c, * knx_c,* kny_c,* knz_c;
  int (*list)[Pm]=new int[Npm][Pm];
  
  x = new double[Np];
  y = new double[Np];
  z = new double[Np];
  x0 = new double[Np];
  y0 = new double[Np];
  z0 = new double[Np];
  x_update = new double[Np];
  y_update = new double[Np];
  z_update = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  vz = new double[Np];
  a = new double[Np];
  kx = new double[Np];
  ky = new double[Np];
  kz = new double[Np];
  knx = new double[Np];
  kny = new double[Np];
  knz = new double[Np];
  kx_w = new double[Np];
  ky_w = new double[Np];
  kz_w = new double[Np];
  knx_w = new double[Np];
  kny_w = new double[Np];
  knz_w = new double[Np];
  kx_c = new double[Np];
  ky_c = new double[Np];
  kz_c = new double[Np];
  knx_c = new double[Np];
  kny_c = new double[Np];
  knz_c = new double[Np];
  kw = new double[Np];
  nx = new double[Np];
  ny = new double[Np];
  nz = new double[Np];
  nx_pr = new double[Np];
  ny_pr = new double[Np];
  nz_pr = new double[Np];
  theta1 = new double[Np];
  theta2 = new double[Np];
  omega = new double[Np];
  r1 = new double[Np];

  char filename[128];
  
  ini_hex(x,y,z,nx,ny,nz,theta1,a,L,Np,kappa);
  initial(nx,ny,nz,Np);
  ini(vx, vy, vz, nx_pr, ny_pr, nz_pr, Np,omega);  
  sprintf(filename,"anc_mix_na_c%.1f_T%.1f.txt",c,temp);
  ofstream file;
  file.open(filename);
  file <<"#"<<" "<<"phi"<<" "<<"rho"<<" "<<"temp"<<" "<<"avU"<<" "<<"cv"<<" "<<"p"<<" "<<"p/T"<<endl;
  
  update(L,Lz,Np,x,y,z,M,Mz,RCHK,list);
  avU0=0.0;
  disp_ave = 0.0;
  
  count=0;
  timer=0.0;
  int count_max=10000;
  int count_coord=0;
  
  copy(x_update,y_update,z_update,x,y,z,Np,x_corr,y_corr,z_corr);
  copy(x0,y0,z0,x,y,z,Np,x_corr,y_corr,z_corr);
  
  for(phi_final = 0.01; phi_final <= 1.0; phi_final+= 0.01){
    avK=0.0,avK0=0.0, p0=0.0,p=0.0,avU=0.0,avU2=0.0;
    for(;;){ // infinite loop
      for(int count=0;count<5000;count++){
	calc_force(x, y, z, L, nx, ny, nz, Np, a, kx, ky, kz, knx,kny,knz,kw,kx_w, ky_w, kz_w, knx_w,kny_w,knz_w,list,theta1,theta2,chi,chi_pr,mu,eta,kappa,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,c,&dummy,&dummy,&dummy,anc);   
	eq_motion(x, y,z, theta1, theta2,vx, vy, vz, nx, ny, nz, nx_pr, ny_pr, nz_pr, omega, dt, kx, ky, kz, knx, kny,knz,kw,kx_w, ky_w, kz_w, knx_w,kny_w,knz_w, Np, &avK0,temp,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,c);
	p_bound(x, y, z, Np, L);
	//cout <<"disp= "<< disp_ave <<" " << "x= "<<" " << x[1]<<" "<< "temp=" << avK0 << "L=" << L <<" "<< "M="<<M << " "<< "phi = " << phi << "count = " << count <<endl;
	
	///////auto update////////////////////
	calc_disp(&disp_max,&disp_ave,x,y,z,x_update,y_update,z_update,x0,y0,z0,Np,L,x_corr,y_corr,z_corr);
	count_num_update++;
	if(disp_max>0.5*0.5){
	  update(L,Lz,Np,x,y,z,M,Mz,RCHK,list);
	  disp_max=0.0;
	  copy(x_update,y_update,z_update,x,y,z,Np,x_corr,y_corr,z_corr);
	  count_num_update=0;
	}
	////////////////////////////
      }
      affine_coord(x, y, z, phi, Np,dphi);  
      phi += dphi;
      L = sqrt(pow(2.,1./2.)*pi*kappa*Np/(6.*kappa*pow(2.,1./6.)*phi));  
      // L = sqrt(pow(2.,1./2.)*pi*kappa*Np/(6.*kappa*phi));
      rho=Np/L/L/lz;
      M=(int)(L/RCHK);       
      if(phi >= phi_final)
        break;
    }
    
    int count2 = 0;
    avK=0.0,avK0=0.0, p0=0.0,p=0.0,avU=0.0,avU2=0.0;
    //count_coord++;
    
    for(;;){
      p0=0.0;
      avU0=0.0;
      calc_force(x, y, z, L, nx, ny, nz, Np, a, kx, ky, kz, knx,kny,knz,kw,kx_w, ky_w, kz_w, knx_w,kny_w,knz_w,list,theta1,theta2,chi,chi_pr,mu,eta,kappa,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,c,&p0,&avU0,&avU2,anc);   
      eq_motion(x, y,z, theta1, theta2,vx, vy, vz, nx, ny, nz, nx_pr, ny_pr, nz_pr, omega, dt, kx, ky, kz, knx, kny,knz,kw,kx_w, ky_w, kz_w, knx_w,kny_w,knz_w, Np, &avK0,temp,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,c);
      p_bound(x, y, z, Np, L);
      avK+=avK0;
      avU+=avU0;
      p += p0;
      //cout <<"disp="<< disp_ave <<"x="<<x[1]<<" "<< "temp=" << avK0 << "L=" << L <<" "<< "M="<<M << " "<< "phi = " << phi << " "<<"count2=" << count2<<endl;
      if(count2 == 100000){
	//if(count_coord==1){ 
	output(x,y,z,nx,ny,nz,x_corr,y_corr,z_corr,L,theta1,theta2,a,Np,temp,c,kappa,phi,anc); 
	 // count_coord=0;}
	file <<phi<<" "<< rho << " " <<avK/count2<<" "<<avU/count2<<" "<<(avU2/count2-avU*avU/count2/count2)/temp/temp<<" "<<avK/count2 + p/count2 <<" " <<(avK/count2 + p/count2)/(avK/count2) <<endl; // to get av
	break;
      }//}
      ///////auto update////////////////////
      calc_disp(&disp_max,&disp_ave,x,y,z,x_update,y_update,z_update,x0,y0,z0,Np,L,x_corr,y_corr,z_corr);
      count2++;
      count_num_update++;
      if(disp_max>0.5*0.5){
	update(L,Lz,Np,x,y,z,M,Mz,RCHK,list);
	disp_max=0.0;
	copy(x_update,y_update,z_update,x,y,z,Np,x_corr,y_corr,z_corr);
	count_num_update=0;
      }
      ////////////////////////////
    }
  }
  file.close();
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] x0;
  delete[] y0;
  delete[] z0;
  delete[] x_update;
  delete[] y_update;
  delete[] z_update;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] a;
  delete[] kx;
  delete[] ky;
  delete[] kz;
  delete[] kx_w;
  delete[] ky_w;
  delete[] kz_w;
  delete[] knx;
  delete[] kny;
  delete[] knz;
  delete[] knx_w;
  delete[] kny_w;
  delete[] knz_w;
  delete[] nx;
  delete[] ny;
  delete[] nz;
  delete[] nx_pr;
  delete[] ny_pr;
  delete[] nz_pr;
  delete[] list;
  delete[] theta1;
  delete[] theta2;
  delete[] omega;
  return 0;
}
