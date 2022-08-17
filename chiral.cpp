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
int copy(double *x_update,double *y_update,double *z_update,double *x,double *y,double *z,int Np,double x_corr,double y_corr,double z_corr){
  int i;
  for(i=0;i<Np;i++){
    x_update[i]=x[i]-x_corr;
    y_update[i]=y[i]-y_corr;
    z_update[i]=z[i]-z_corr;
  }
  return 0;
}

int calc_disp(double *disp_max,double *disp_ave,double *x,double *y,double *z,double *x_update,double *y_update,double *z_update,double *x0,double *y0,double *z0,int Np,int L,double x_corr,double y_corr,double z_corr){
  int i;
  double dx,dy,dz;
  double disp;
  
  *disp_ave = 0.0;

  for(i=0;i<Np;i++){
    dx=x[i]-x_update[i]-x_corr;
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    dy=y[i]-y_update[i]-y_corr;
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    dz=z[i]-z_update[i]-z_corr;
    if(dz>0.5*L) dz-=L;
    else if(dz<-0.5*L)dz+=L;
    disp = dx*dx+dy*dy+dz*dz;
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
    dz=z[i]-z0[i]-z_corr;
    if(dz>0.5*L) dz-=L;
    else if(dz<-0.5*L)dz+=L;
    disp = dx*dx+dy*dy+dz*dz;
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
    if(dz>0.5*L) dz-=L;
    else if(dz<-0.5*L)dz+=L;

    disp = dx*dx+dy*dy+dz*dz;    
    if(disp > *disp_max)
      *disp_max =disp;    
  }
  
  return 0;
}


int com_correction(double *x,double *y,double *z,double *x_corr,double *y_corr,double *z_corr,int Np,double L){
  int i;
  double dx,dy,dz;
  static double x0[Npm],y0[Npm],z0[Npm];
  static bool IsFirst = true;
  if(IsFirst){
    for(i=0;i<Np;i++){
      x0[i]=x[i];
      y0[i]=y[i];
      z0[i]=z[i];
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

    dz=z[i]-z0[i];
    if(dz>0.5*L) dz-=L;
    else if(dz<-0.5*L)dz+=L;
    
    *x_corr+=dx/Np; //center of mass displacement.x
    *y_corr+=dy/Np;
    *z_corr+=dz/Np;
    
    x0[i]=x[i];
    y0[i]=y[i];
    z0[i]=z[i];
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
    if (z[k] < 0.0) {
      z[k] = z[k] + L;
    }    
    if (z[k] >  L) {
      z[k] = z[k] - L;
    }
  }
  return 0;
}

int output(double *x,double *y,double *z,double* nx, double* ny,double* nz,double x_corr,double y_corr,double z_corr,double L,double *theta1,double *theta2,double *r1,int Np){
  int i;
  static int count_file=1;
  double x1[Npm],y1[Npm],z1[Npm];
  char filename[128];   
  sprintf(filename,"coord_0.6_1.0_rand_affine.dat");
  ofstream file;
  file.open(filename);
  for(i=0;i<Np;i++){
    x1[i]=x[i]-x_corr;
    y1[i]=y[i]-y_corr;
    z1[i]=z[i]-z_corr;
    //theta1[i] = acos(nz[i]);
    //theta2[i] = atan(ny[i]/(nx[i]+0.000000001));

  }
  p_bound(x1, y1, z1, Np, L);  
  for(i=0;i<Np;i++)
    file << x1[i] << " " << y1[i] << " " << z1[i]<< " "<<nx[i]<<" " <<ny[i]<<" " << nz[i] << endl; 
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
void affine_coord(double* x,double* y,double* z,double phi,int Np){
  for(int k = 0; k < Np; k++){
    x[k] = pow(phi/(phi + 0.001),1./3)*x[k];
    y[k] = pow(phi/(phi + 0.001),1./3)*y[k];
    z[k] = pow(phi/(phi + 0.001),1./3)*z[k];
  }
}




void ini_fcc(double *x,double *y,double*z,double* nx,double* ny, double* nz, double* theta1, double* theta2, double *a,double L,int Np){
  int num_x = (int)(pow(Np/4.,1./3.))+1;
  int num_y = (int)(pow(Np/4.,1./3.))+1;
  int num_z = (int)(pow(Np/4.,1./3.))+1;
  double c[3*Np+3];
  for(int j =0 ; j< 3*Np +3; j++){
    c[j] = 0.; 
  }
  

  double lx = L/double(num_x);
  double ly = L/double(num_y);
  double lz = L/double(num_z);

 // cout << "lx" <<" " <<lx << endl;

  int i,j,k=0;
  int p = 0;
  for(k=0;k<num_z;k++){
    for(j=0;j<num_y;j++){
      for(i=0;i<num_x;i++){
        c[p] = i*lx;
        c[p+1] = j*ly;
        c[p+2] = k*lz;
    //    cout <<  "p = " << "  " << p << " " << ", x = " << "  " << c[p]<< " " << ", y = " << "  " << c[p+1] << " " << ", z = " << "  " << c[p+2] << endl;
        p += 3;
     
        if(p >=3*Np+3)
        break;

        c[p] = i*lx;
        c[p+1] = j*ly + ly/sqrt(2.);
        c[p+2] = k*lz + lz/sqrt(2.);
    //    cout <<  "p = " << "  " << p << " " << ", x = " << "  " << c[p]<< " " << ", y = " << "  " << c[p+1] << " " << ", z = " << "  " << c[p+2] << endl;
        p += 3;
   
        if(p >=3*Np+3)
        break;

        c[p] = i*lx + lz/sqrt(2.);
        c[p+1] = j*ly;
        c[p+2] = k*lz + lz/sqrt(2.);
    //    cout <<  "p = " << "  " << p << " " << ", x = " << "  " << c[p]<< " " << ", y = " << "  " << c[p+1] << " " << ", z = " << "  " << c[p+2] << endl;
        p += 3;
        if(p >=3*Np+3)
        break;

        c[p] = i*lx + lz/sqrt(2.);
        c[p+1] = j*ly + ly/sqrt(2.);
        c[p+2] = k*lz;
    //    cout <<  "p = " << "  " << p << " " << ", x = " << "  " << c[p]<< " " << ", y = " << "  " << c[p+1] << " " << ", z = " << "  " << c[p+2] << endl;
        p += 3;
        if(p >=3*Np+3)
        break;
       }
        if(p >=3*Np+3)
        break;       
        
      }
        if(p >=3*Np+3)
        break;
    }

    for(p = 0 ; p < Np; p++){
        x[p] = c[3*p];
        y[p] = c[3*p+1];
        z[p] = c[3*p+2];
        nx[p] = unif_rand(-1, 1);
        ny[p] = unif_rand(-1, 1);
        nz[p] = unif_rand(-1, 1);
        a[p] = 1.;
       // theta1[p] = 1./sqrt(2.);
       // theta2[p] = 1./sqrt(2.);
    }
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

void update(double L,int Np,double *x,double *y,double *z,int M,double Rcell,int (* list)[Pm]) 
{
  int i,j,k;
  int nx,ny,nz;
  int l,m,n;
  double dx,dy,dz,r;
  
  int (*map)[Npm]=new int[M*M*M][Npm];
  for(k=0;k<M;k++)
    for(j=0;j<M;j++)
     for(i=0;i<M;i++)
     	 map[i+M*j+M*M*k][0]=0;
  for(i=0;i<Np;i++){
   nx=f((int)(x[i]*M/L),M);
   ny=f((int)(y[i]*M/L),M);
   nz=f((int)(z[i]*M/L),M); //address of particle i is determined. 
   for(n=nz-1;n<=nz+1;n++){
    for(m=ny-1;m<=ny+1;m++){
     for(l=nx-1;l<=nx+1;l++){
       map[f(l,M)+M*f(m,M)+M*M*f(n,M)][map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0] +1]=i; 
       map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0]++;
        }
       } 
     }
}
  
for (i=0;i<Np;i++){ 
  list[i][0]=0;
  nx = f((int)(x[i]*M/L),M); 
  ny = f((int)(y[i]*M/L),M); 
  nz = f((int)(z[i]*M/L),M);
    
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
      if(dz<-L/2.0)
       dz+=L;
      else if(dz> L/2.0)
       dz-=L;

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

int calc_force(double* x, double* y, double* z, double L, double* nx, double* ny, double* nz,int Np, double* a, double* kx, double* ky,double* kz,double *knx,double *kny,double *knz,double* kx_c, double* ky_c,double* kz_c,double *knx_c,double *kny_c,double *knz_c,int (*list)[Pm],double *theta1,double *theta2,double chi,double chi_pr,double mu,double eta) {
  int i, j, k;
  // *avU = 0.0;
  double r,r2,r3,r4,ri,rj,cij,dij;
  double R2p,R2n,R1p,R1n,R2p_pr,R2n_pr,R1p_pr,R1n_pr;
  double t,drU,diU,djU,dcU,U,drUc,diUc,djUc,dcUc,ddUc,Uc,A,B,A_c,B_c,A1,A2,A6,A12,fx_ij,fy_ij,fz_ij,fx_ij_c,fy_ij_c,fz_ij_c,U1;
  double e1,e2;
  double dx, dy, dz;
  double aij, aij_pr, aij3, aij_pr3;
  double cut;
  // cut = 2.0;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    kz[k] = 0.0;
    knx[k] = 0.0;
    kny[k] = 0.0;
    knz[k] = 0.0;
    kx_c[k] = 0.0;
    ky_c[k] = 0.0;
    kz_c[k] = 0.0;
    knx_c[k] = 0.0;
    kny_c[k] = 0.0;
    knz_c[k] = 0.0;
  }
 for (i = 0; i < Np; i++)
    {
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
	  if(r2< 3.0*3.0){  // the force by potential energy would work if the distance is smaller than 3*(charcteristic length) maybe i better put it just after calculating r??
	    
	    ri = dx*nx[i] + dy*ny[i] + dz*nz[i];
	    rj = dx*nx[list[i][j]] + dy*ny[list[i][j]] + dz*nz[list[i][j]];
	    cij = nx[i]*nx[list[i][j]] + ny[i]*ny[list[i][j]] + nz[i]*nz[list[i][j]];
        dij = dx*(ny[i]*nz[list[i][j]] - nz[i]*ny[list[i][j]]) + dy*(nz[i]*nx[list[i][j]] - nx[i]*nz[list[i][j]]) + dz*(nx[i]*ny[list[i][j]] - ny[i]*nx[list[i][j]]);
	    
	    
	    R2p =(ri+rj)*(ri+rj)/(1+chi*cij);
 	    R2n =(ri-rj)*(ri-rj)/(1-chi*cij);
 	    R1p =(ri+rj)/(1+chi*cij);
 	    R1n =(ri-rj)/(1-chi*cij);

 	    R2p_pr =(ri+rj)*(ri+rj)/(1+chi_pr*cij);
 	    R2n_pr =(ri-rj)*(ri-rj)/(1-chi_pr*cij);
 	    R1p_pr =(ri+rj)/(1+chi_pr*cij);
 	    R1n_pr =(ri-rj)/(1-chi_pr*cij);    
	    	    
	    aij_pr = aij/sqrt(1.-(chi/2./r2)*(R2p+R2n));
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
	    
        // Gay Berne //
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

	    drUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*cij*(((A_c*mu*chi_pr)/e2/r4)*(R2p_pr + R2n_pr) - ((aij_pr3*chi*B_c)/2./aij3/r4)*(R2p + R2n) - B_c/(aij*r) - A_c/r2);//analytical calculation of the 1'st derivative / r
	    diUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*cij*(-((A_c*mu*chi_pr)/e2/r3)*(R1p_pr + R1n_pr)+((aij_pr3*chi*B_c)/2./aij3/r3)*(R1p + R1n));//analytical calculation of the 1'st derivative / r*ni
        djUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*cij*(-((A_c*mu*chi_pr)/e2/r3)*(R1p_pr - R1n_pr)+((aij_pr3*chi*B_c)/2./aij3/r3)*(R1p - R1n));//analytical calculation of the 1'st derivative / r*ni
        dcUc = 4.*pow(e1,eta)*pow(e2,mu)*dij*(((A_c*mu*chi_pr*chi_pr)/2./e2/r3)*(R1p_pr*R1p_pr - R1n_pr*R1n_pr) - ((aij_pr3*chi*chi*B_c)/4./aij3/r3)*(R1p*R1p - R1n*R1n) + A_c*eta*chi*chi*cij*cij*e1*e1/r);//analytical calculation of the 1'st derivative / ri*ni
        ddUc = 4.*pow(e1,eta)*pow(e2,mu)*cij*A_c/r;

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
        kny_c[list[i][j]] -= djUc*(dy - rj*ny[list[i][j]]) + dcUc*(ny[i] - cij*ny[list[i][j]]) + ddUc*(-ri*(-nz[list[i][j]]*nx[i] + nx[list[i][j]]*nz[i]) + cij*(dz*nx[list[i][j]] - dx*nz[list[i][j]]));

	    knz_c[i] -= diUc*(dz - ri*nz[i]) + dcUc*(nz[list[i][j]] - cij*nz[i]) + ddUc*(ri*(nx[list[i][j]]*ny[i] - ny[list[i][j]]*nx[i]) - cij*(dx*ny[i] - dy*nx[i]));
        knz_c[list[i][j]] -= djUc*(dz - rj*nz[list[i][j]]) + dcUc*(nz[i] - cij*nz[list[i][j]]) + ddUc*(-ri*(-nx[list[i][j]]*ny[i] + ny[list[i][j]]*nx[i]) + cij*(dx*ny[list[i][j]] - dy*nx[list[i][j]]));


	    U1=4.*pow(e1,eta)*pow(e2,mu)*A + 4.*pow(e1,eta)*pow(e2,mu)*cij*dij*A_c/r;
	  }
        }
    }
 
 return 0; 
}



int eq_motion(double* x, double* y,double* z,double* theta1,double* theta2, double* vx, double* vy, double* vz,double* nx,double* ny,double* nz, double* nx_pr,double* ny_pr,double* nz_pr, double* omega, double dt, double* kx, double* ky, double* kz,double* knx,double* kny,double* knz,double* kx_c, double* ky_c, double* kz_c,double* knx_c,double* kny_c,double* knz_c, double c, int Np, double* avK, double Th) {
  double zeta, zeta_r;
  double n[Np];
  double n_pr2[Np];
  double Rx[Np];
  double Ry[Np];
  double Rz[Np];

  int k;
for (k = 0; k < Np; k++) {
    n[k] = 0.0;
    n_pr2[k] = 0.0;
    Rx[k] = 0.0;
    Ry[k] = 0.0;
    Rz[k] = 0.0;
  }

  zeta = 1.;
  zeta_r = 100.;
  *avK=0.0;
  for (k = 0; k < Np; k++) {
    vx[k] += -vx[k] * zeta * dt + kx[k] * dt + c * kx_c[k] + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vy[k] += -vy[k] * zeta * dt + ky[k] * dt + c * ky_c[k] + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vz[k] += -vz[k] * zeta * dt + kz[k] * dt + c * kz_c[k] + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    Rx[k] = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    Ry[k] = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    Rz[k] = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    n_pr2[k] = nx_pr[k]*nx_pr[k] + ny_pr[k]*ny_pr[k] + nz_pr[k]*nz_pr[k];
    nx_pr[k] += -n_pr2[k]*nx[k]*dt - zeta_r*nx_pr[k]*dt + knx[k]*dt + c*knx[k]*dt + Ry[k]*nz[k] - Rz[k]*ny[k] ;
    ny_pr[k] += -n_pr2[k]*ny[k]*dt - zeta_r*ny_pr[k]*dt + kny[k]*dt + c*kny[k]*dt + Rz[k]*nx[k] - Rx[k]*nz[k] ;
    nz_pr[k] += -n_pr2[k]*nz[k]*dt - zeta_r*nz_pr[k]*dt + knz[k]*dt + c*knz[k]*dt + Rx[k]*ny[k] - Ry[k]*nx[k] ;
  
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    z[k] += vz[k] * dt;

    nx[k] += nx_pr[k] * dt;
    ny[k] += ny_pr[k] * dt;
    nz[k] += nz_pr[k] * dt;

    n[k] = sqrt(nx[k] * nx[k] + ny[k] * ny[k] + nz[k] * nz[k]);

    nx[k] = nx[k] / n[k];
    ny[k] = ny[k] / n[k];
    nz[k] = nz[k] / n[k];

//    theta1[k] = acos(nz[k]);
 //   theta2[k] = acos(nx[k]/(sin(theta1[k])+0.000000001));

   
    *avK += vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k] + (nx_pr[k] * nx_pr[k] + ny_pr[k] * ny_pr[k] + nz_pr[k] * nz_pr[k] )/(nx[k] * nx[k] + ny[k] * ny[k] + nz[k] * nz[k]);
  }
  *avK = *avK / Np / 2.0/ 2.5;
  return 0;
}

int main(int argc, char *argv[])
{
  srand((unsigned) time(NULL));
  double t,avU=0.0,avU2=0.0,avU0=0.0,avK=0.0,avK0=0.0,dummy;
  double x_corr=0.0,y_corr=0.0,z_corr=0.0;
  int i,count=0,count_num_update=0;

  double disp_max=0.0;
  double disp_ave;
  int Np = 1000;int count_th=200;
  double disp_th2 = 180;
  double dt =0.001;//  //parameters;
  double temp = 0.6;
  double Th;
  double phi=0.01;


  double kappa = 3.; // shape anisotropy parameter
  double kappa_pr = 3.; // energy anisotropy parameter
  double chi = (kappa*kappa-1.)/(kappa*kappa+1.); //atof(argv[1]);  
 
  double eta = 1.0;
  double mu = 2.0;
  double chi_pr = (pow(kappa_pr,1./mu)-1.)/(pow(kappa_pr,1./mu)+1.) ; 

  double timer;

  double RCHK=4.0;// cut off +1.0
  double L;
         L = pow((pi*pow(2.,1./2.)*kappa*Np)/(6.*phi),1./3.);//(r1+r2)/4.0*sqrt(pow(2.,1./3.)*pi*kappa*Np/phi);
  int M;
         M=(int)(L/RCHK);
  cout << "L=" << L <<" "<< "M="<<M <<endl;
  
  double* x, * y, *z,* x0, * y0, *z0,* vx, * vy,*vz, * a, * kx, * ky,*kz, * knx,* kny,* knz,* kx_c, * ky_c,*kz_c, * knx_c,* kny_c,* knz_c, * nx,* ny,* nz, * nx_pr,* ny_pr,* nz_pr,*x_update,*y_update,*z_update, *theta1, *theta2 ,*omega;
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
  kx_c = new double[Np];
  ky_c = new double[Np];
  kz_c = new double[Np];
  knx_c = new double[Np];
  kny_c = new double[Np];
  knz_c = new double[Np];
  nx = new double[Np];
  ny = new double[Np];
  nz = new double[Np];
  nx_pr = new double[Np];
  ny_pr = new double[Np];
  nz_pr = new double[Np];
  theta1 = new double[Np];
  theta2 = new double[Np];
  omega = new double[Np];

  char filename[128];
  
  // ini_coord_square(x,y,theta,a,L,Np);
  //ini_hex(x,y,theta,a,L,Np);
 
  ini_fcc(x,y,z,nx,ny,nz,theta1,theta2,a,L,Np);

  ini(vx, vy, vz, nx_pr, ny_pr, nz_pr, Np,omega);  

  
  sprintf(filename,"%s/energy_time.txt",argv[1]);
  ofstream file;
  file.open(filename);
  
// cout << "before cell list" << endl;
  update(L,Np,x,y,z,M,RCHK,list);
  avU0=0.0;
  disp_ave = 0.0;
// cout << "done cell list" << endl;

  count=0;
  timer=0.0;
  
  copy(x_update,y_update,z_update,x,y,z,Np,x_corr,y_corr,z_corr);
  copy(x0,y0,z0,x,y,z,Np,x_corr,y_corr,z_corr);

double c = 1.0 ;// chirality
for(;;){ // infinite loop

   for(int count=0;count<1000;count++){
    calc_force(x, y, z, L, nx, ny, nz, Np, a, kx, ky, kz, knx,kny,knz,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,list,theta1,theta2,chi,chi_pr,mu,eta);   
    eq_motion(x, y,z, theta1, theta2,vx, vy, vz, nx, ny, nz, nx_pr, ny_pr, nz_pr, omega, dt, kx, ky, kz, knx, kny,knz,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,c,Np, &avK0,temp);
    com_correction(x,y,z,&x_corr,&y_corr,&z_corr, Np, L);
    p_bound(x, y, z, Np, L);
    cout <<"disp= "<< disp_ave <<" " << "x= "<<" " << x[1]<<" "<< "temp=" << avK0 << "L=" << L <<" "<< "M="<<M << " "<< "phi = " << phi << "count = " << count <<endl;
    count++;
       ///////auto update////////////////////
    calc_disp(&disp_max,&disp_ave,x,y,z,x_update,y_update,z_update,x0,y0,z0,Np,L,x_corr,y_corr,z_corr);
    count_num_update++;
    if(disp_max>0.5*0.5){
      update(L,Np,x,y,z,M,RCHK,list);
      disp_max=0.0;
      copy(x_update,y_update,z_update,x,y,z,Np,x_corr,y_corr,z_corr);
      count_num_update=0;
    }
    ////////////////////////////
    }
    
    affine_coord(x, y, z, phi, Np);  
       phi += 0.001;
     //phi += 0.001*exp(-phi);   
      L = pow((pi*pow(2.,1./2.)*kappa*Np)/(6.*phi),1./3.);
      M=(int)(L/RCHK); 
    
    if(phi>= 1.0)
    break;}


    disp_ave = 0.0;
    int count2 = 0;
    for(;;){
    calc_force(x, y, z, L, nx, ny, nz, Np, a, kx, ky, kz, knx,kny,knz,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,list,theta1,theta2,chi,chi_pr,mu,eta);   
    eq_motion(x, y,z, theta1, theta2,vx, vy, vz, nx, ny, nz, nx_pr, ny_pr, nz_pr, omega, dt, kx, ky, kz, knx, kny,knz,kx_c, ky_c, kz_c, knx_c,kny_c,knz_c,c,Np, &avK0,temp);
    com_correction(x,y,z,&x_corr,&y_corr,&z_corr, Np, L);
    p_bound(x, y, z, Np, L);
    cout <<"disp="<< disp_ave <<"x="<<x[1]<<" "<< "temp=" << avK0 << "L=" << L <<" "<< "M="<<M << " "<< "phi = " << phi << " "<<"count2=" << count2<<endl;
     if(count2 > 100000){
    // if(disp_ave>disp_th2){ 
      output(x,y,z,nx,ny,nz,x_corr,y_corr,z_corr,L,theta1,theta2,a,Np); 
      break;
      }
   ///////auto update////////////////////
    calc_disp(&disp_max,&disp_ave,x,y,z,x_update,y_update,z_update,x0,y0,z0,Np,L,x_corr,y_corr,z_corr);
     count2++;
    count_num_update++;
    if(disp_max>0.5*0.5){
      update(L,Np,x,y,z,M,RCHK,list);
      disp_max=0.0;
      copy(x_update,y_update,z_update,x,y,z,Np,x_corr,y_corr,z_corr);
      count_num_update=0;
    }
    ////////////////////////////
    
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