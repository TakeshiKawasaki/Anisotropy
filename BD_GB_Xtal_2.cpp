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

int output(double *x,double *y,double x_corr,double y_corr,double L,double *theta,double *r1,int Np){
  int i;
  static int count_file=1;
  double x1[Npm],y1[Npm];
  char filename[128];   
  sprintf(filename,"coord_%.d.dat",count_file);
  ofstream file;
  file.open(filename);
  for(i=0;i<Np;i++){
    x1[i]=x[i]-x_corr;
    y1[i]=y[i]-y_corr;
  }
  p_bound(x1, y1, Np, L);  
  for(i=0;i<Np;i++)
    file << x1[i] << " " << y1[i]<< " "<<theta[i]<<" " << r1[i] << endl; 
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


void ini_coord_square(double *x,double *y,double *theta,double *a,double L,int Np){
  int num_x = (int)sqrt(Np)+1;
  int num_y = (int)sqrt(Np)+1;
  int i,j,k=0;
  for(j=0;j<num_y;j++){
    for(i=0;i<num_x;i++){
      x[i+num_x*j] = i*L/(double)num_x;
      y[i+num_x*j] = j*L/(double)num_y;
      theta[i+num_x*j]=M_PI/4.0;
      a[i+num_x*j]=1.0;
      k++;
      if(k>=Np)
        break;
    }
    if(k>=Np)
      break;
  }
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

int calc_force(double* x, double* y, double L, int Np, double* a, double* kx, double* ky,double *kth,int (*list)[Pm],double *theta,double chi,double chi_pr,double mu,double eta) {
  int i, j, k;
  // *avU = 0.0;
  double r,r2,r4,ri,rj,cij;
  double R2p,R2n,R1p,R1n,R2p_pr,R2n_pr,R1p_pr,R1n_pr;
  double t,drU,diU,djU,dcU,U,A,B,A1,A2,A6,A12,fx_ij,fy_ij,U1;
  double e1,e2;
  double dx, dy;
  double aij, aij_pr;
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
	  // cout<<theta[list[i][j]]<<endl;
          r = sqrt(dx*dx+dy*dy);
          r2 = dx*dx+dy*dy;
          r4 = r2*r2;
	  if(r2< 7.0*7.0){  // the force by potential energy would work if the distance is smaller than 3*(charcteristic length) maybe i better put it just after calculating r??
	    
	    ri = dx*cos(theta[i])+dy*sin(theta[i]);
	    rj = dx*cos(theta[list[i][j]])+dy*sin(theta[list[i][j]]);
	    cij = cos(theta[i])*cos(theta[list[i][j]])+sin(theta[i])*sin(theta[list[i][j]]);
	    
	    
	    R2p =(ri+rj)*(ri+rj)/(1+chi*cij);
 	    R2n =(ri-rj)*(ri-rj)/(1-chi*cij);
 	    R1p =(ri+rj)/(1+chi*cij);
 	    R1n =(ri-rj)/(1-chi*cij);

 	    R2p_pr =(ri+rj)*(ri+rj)/(1+chi_pr*cij);
 	    R2n_pr =(ri-rj)*(ri-rj)/(1-chi_pr*cij);
 	    R1p_pr =(ri+rj)/(1+chi_pr*cij);
 	    R1n_pr =(ri-rj)/(1-chi_pr*cij);    
	    	    
	    aij_pr = aij/sqrt(1.-(chi/2./r2)*(R2p+R2n));
	    
	    A1 = (aij/(r-aij_pr+aij));
	    A2 =A1*A1;
	    A6 = A2*A2*A2;
	    A12 = A6*A6;
	    
	    A = A12-A6;
	    B = 13.*A12*A1-7.*A6*A1;
	    
	    
	    //	    e1 = 1/sqrt(1-chi*chi-cij*cij);  //epsilon1
	    e1 = 1./sqrt(1-chi*chi*cij*cij);  //epsilon1
	    e2 = 1.-(chi_pr/2./r2)*(R2p_pr+R2n_pr); // epsilon2
	    
	    
	    //start derivation
	    
	    drU = 4.*pow(e1,eta)*pow(e2,mu)*(((A*mu*chi_pr)/e2/r4)*(R2p_pr+R2n_pr)-((aij_pr*aij_pr*aij_pr*chi*B)/2./r4)*(R2p+R2n)-B/r);//analytical calculation of the 1'st derivative / r
	    diU = 4.*pow(e1,eta)*pow(e2,mu)*(-((A*mu*chi_pr)/e2/r2)*(R1p_pr+R1n_pr)-((aij_pr*aij_pr*aij_pr*chi*B)/2./r2)*(R1p+R1n));//analytical calculation of the 1'st derivative / r*ni
        djU = 4.*pow(e1,eta)*pow(e2,mu)*(-((A*mu*chi_pr)/e2/r2)*(R1p_pr-R1n_pr)-((aij_pr*aij_pr*aij_pr*chi*B)/2./r2)*(R1p-R1n));//analytical calculation of the 1'st derivative / r*nj
        dcU = 4.*pow(e1,eta)*pow(e2,mu)*(((A*mu*chi_pr*chi_pr)/2./e2/r2)*(R1p_pr*R1p_pr-R1n_pr*R1n_pr)-((aij_pr*aij_pr*aij_pr*chi*chi*B)/4./r2)*(R1p*R1p-R1n*R1p)+A*eta*chi*chi*cij*e1*e1);//analytical calculation of the 1'st derivative / ri*ni
    
	    fx_ij= drU*dx/r+diU*cos(theta[i])+djU*cos(theta[list[i][j]]);
	    fy_ij= drU*dy/r+diU*sin(theta[i])+djU*sin(theta[list[i][j]]);
	    
	    kx[i]+=fx_ij;
	    kx[list[i][j]]-=fx_ij;
	    
	    ky[i]+=fy_ij;	    
	    ky[list[i][j]]-=fy_ij;
	    
	    kth[i]         -=diU*(cos(theta[i])*dy-sin(theta[i])*dx)+dcU*(cos(theta[i])*sin(theta[list[i][j]])-sin(theta[i])*cos(theta[list[i][j]]));  
	    kth[list[i][j]]-=djU*(cos(theta[list[i][j]])*dy-sin(theta[list[i][j]])*dx)+dcU*(cos(theta[list[i][j]])*sin(theta[i])-sin(theta[list[i][j]])*cos(theta[i]));
	    U1=4.*e1*e2*A;
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
  double t,avU=0.0,avU2=0.0,avU0=0.0,avK=0.0,avK0=0.0,dummy;
  double x_corr=0.0,y_corr=0.0;
  int i,count=0,count_num_update=0;

  double disp_max=0.0;
  double disp_ave=0.0;
  int Np = 1000;int count_th=200;
  double r1=1.0, r2=1.0;
  double disp_th =10.0;
  double dt =0.0005;//  //parameters;
  double time_stable_1 = 100.;
  double time_stable_2 = 1000.;
  double Th;


  double kappa=2.; // shape anisotropy parameter
  double kappa_pr=2.; // energy anisotropy parameter
  double chi = (kappa*kappa-1.)/(kappa*kappa+1.); //atof(argv[1]);  
 
  double eta = 1.0;
  double mu = 2.0;
  double chi_pr = (pow(kappa_pr*kappa_pr,1./mu)-1.)/(pow(kappa_pr*kappa_pr,1./mu)+1.) ; 

  double timer;
  double chi0=0.2;

  double RCHK=9.0;
  double L = (r1+r2)/4.0*sqrt(pi*kappa*Np/0.35);
  int    M=(int)(L/RCHK);
  cout << "L=" << L <<" "<< "M="<<M <<endl;
  
  double* x, * y,* x0, * y0, * vx, * vy, * a, * kx, * ky,*kth,*x_update,*y_update,*theta,*omega;
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
  theta = new double[Np];
  omega = new double[Np];

  char filename[128];
  
  ini_coord_square(x,y,theta,a,L,Np);
  ini(vx, vy, Np,omega);  
  
  sprintf(filename,"%s/energy_time.txt",argv[1]);
  ofstream file;
  file.open(filename);
  
  update(L,Np,x,y,M,RCHK,list);
  avU0=0.0;


  count=0;
  timer=0.0;
  
  copy(x_update,y_update,x,y,Np,x_corr,y_corr);
  copy(x0,y0,x,y,Np,x_corr,y_corr);
  
  Th=0.001;
  for(;;){ // infinite loop
    calc_force(x, y, L, Np, a, kx, ky,kth,list,theta,chi,chi_pr,mu,eta);
    eq_motion(x, y, theta, vx, vy, omega, dt, kx, ky,kth, Np, &avK0,Th);
    com_correction(x,y,&x_corr,&y_corr,Np, L);
    p_bound(x, y, Np, L);
    
    cout <<"x="<<x[1]<<" "<<avK0<<endl;      
    if(disp_ave>disp_th*disp_th){ // the condition for changing temp
      output(x,y,x_corr,y_corr,L,theta,a,Np); // does it mean every output has different temp??
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
   
  delete[] kth;
  delete[] list;
  delete[] theta;
  delete[] omega;
  return 0;
}
