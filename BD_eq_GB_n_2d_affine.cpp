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

int output(double *x,double *y, double* nx,double* ny,double x_corr,double y_corr,double L,double *theta,double *r1,int Np){
  int i;
  static int count_file=1;
  double x1[Npm],y1[Npm];
  char filename[128];   
  sprintf(filename,"coord_1.0_1.0_10000_rand_affine.dat");
  ofstream file;
  file.open(filename);
  for(i=0;i<Np;i++){
    x1[i]=x[i]-x_corr;
    y1[i]=y[i]-y_corr;
    theta[i] = atan(ny[i]/(nx[i]+0.000001));
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

void ini_hex(double *x,double *y,double* nx, double* ny,double *theta,double *a,double L,int Np){
  int num_x = (int)sqrt(Np)+1;
  int num_y = (int)sqrt(Np)+1;
  int i,j,k=0;
  double shift;
  for(j=0;j<num_y;j++){
    for(i=0;i<num_x;i++){
      shift=(double)j*0.5-j/2;
      x[i+num_x*j] = (shift+i)*L/(double)num_x;
      y[i+num_x*j] = j*L/(double)num_y;
      nx[i+num_x*j] = unif_rand(-1, 1);
      ny[i+num_x*j] = unif_rand(-1, 1);
//      theta[i+num_x*j] = atan(ny[i+num_x*j]/nx[i+num_x*j]);
      a[i+num_x*j]=1.0;
      k++;
      if(k>=Np)
        break;
    }
    if(k>=Np)
      break;
  }
//  char filename[128];   
//  sprintf(filename,"coord_1_0.3_0.3_rand_ini.dat");
//  ofstream file;
//  file.open(filename);
//  for(k = 0; k <Np; k++){
//    file << x[k] << " " <<y[k]<< " "<<theta[k]<<" " << a[k] << endl; 
//  }
//  file.close();
}  

void affine_coord(double *x,double *y, double phi, int Np, double dphi) {
  for(int k = 0; k < Np; k++){
    //x[k] = sqrt(phi/(phi + 0.001*exp(-phi)))*x[k];
   // y[k] = sqrt(phi/(phi + 0.001*exp(-phi)))*y[k];
    x[k] = sqrt(phi/(phi + dphi))*x[k];
    y[k] = sqrt(phi/(phi + dphi))*y[k];
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

int ini(double* vx, double* vy, double* nx_pr, double* ny_pr, int Np, double* omega) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
    nx_pr[j] = 0.0;
    ny_pr[j] = 0.0;
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

int calc_force(double* x, double* y, double* nx, double* ny, double L, int Np, double* a, double* kx, double* ky,double *knx, double* kny, double* kth,int (*list)[Pm],double *theta,double chi,double chi_pr,double mu,double eta) {
  int i, j, k;
  // *avU = 0.0;
  double r,r2,r3,r4,ri,rj,cij;
  double R2p,R2n,R1p,R1n,R2p_pr,R2n_pr,R1p_pr,R1n_pr;
  double t,drU,diU,djU,dcU,U,A,B,A1,A2,A6,A12,fx_ij,fy_ij,U1;
  double e1,e2;
  double dx, dy;
  double aij, aij_pr, aij3, aij_pr3;
  double cut;
  // cut = 2.0;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    knx[k] = 0.0;
    kny[k] = 0.0;
    kth[k] = 0.0;
  }
 for (i = 0; i < Np; i++)
    {
      for (j = 1; j <=list[i][0]; j++)
	{
          dx = x[i] - x[list[i][j]];
          dy = y[i] - y[list[i][j]];
          if (dx > (0.5 * L))
            dx -= L;
          if (dx < -(0.5 * L))
            dx += L;
          if (dy > (0.5 * L))
            dy -= L;
          if (dy < -(0.5 * L))
            dy += L;

          aij = (a[list[i][j]] + a[i]) / 2.0;
          aij3 = aij*aij*aij;

	  // cout<<theta[list[i][j]]<<endl;
          r = sqrt(dx*dx+dy*dy);
          r2 = dx*dx+dy*dy;
          r3 = r2*r;
          r4 = r2*r2;

	  if(r2< 3.0*3.0){  // the force by potential energy would work if the distance is smaller than 3*(charcteristic length) maybe i better put it just after calculating r??
	    
	    ri = dx*nx[i]+dy*ny[i];
	    rj = dx*nx[list[i][j]]+dy*ny[list[i][j]];
	    cij = nx[i]*nx[list[i][j]]+ny[i]*ny[list[i][j]];
	    
	    
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
	    
	    
	    //	    e1 = 1/sqrt(1-chi*chi-cij*cij);  //epsilon1
	    e1 = 1./sqrt(1-chi*chi*cij*cij);  //epsilon1
	    e2 = 1.-(chi_pr/2./r2)*(R2p_pr+R2n_pr); // epsilon2
	    
	    
	    //start derivation
	    
	    drU = 4.*pow(e1,eta)*pow(e2,mu)*(((A*mu*chi_pr)/e2/r3)*(R2p_pr+R2n_pr)-((aij_pr3*chi*B)/2./aij3/r3)*(R2p+R2n)-B/aij);//analytical calculation of the 1'st derivative / r
	    diU = 4.*pow(e1,eta)*pow(e2,mu)*(-((A*mu*chi_pr)/e2/r2)*(R1p_pr+R1n_pr)+((aij_pr3*chi*B)/2./aij3/r2)*(R1p+R1n));//analytical calculation of the 1'st derivative / r*ni
        djU = 4.*pow(e1,eta)*pow(e2,mu)*(-((A*mu*chi_pr)/e2/r2)*(R1p_pr-R1n_pr)+((aij_pr3*chi*B)/2./aij3/r2)*(R1p-R1n));//analytical calculation of the 1'st derivative / r*nj
        dcU = 4.*pow(e1,eta)*pow(e2,mu)*(((A*mu*chi_pr*chi_pr)/2./e2/r2)*(R1p_pr*R1p_pr-R1n_pr*R1n_pr)-((aij_pr3*chi*chi*B)/4./aij3/r2)*(R1p*R1p-R1n*R1n)+A*eta*chi*chi*cij*e1*e1);//analytical calculation of the 1'st derivative / ri*ni
    
	    fx_ij= drU*dx/r+diU*nx[i]+djU*nx[list[i][j]];
	    fy_ij= drU*dy/r+diU*ny[i]+djU*ny[list[i][j]];
	    
	    kx[i]-=fx_ij;
	    kx[list[i][j]]+=fx_ij;
	    
	    ky[i]-=fy_ij;	    
	    ky[list[i][j]]+=fy_ij;
	    
//      knx[i] -= diU*dx + dcU*nx[list[i][j]];
//      knx[list[i][j]] -= djU*dx + dcU*nx[i];

//    kny[i] -= diU*dy + dcU*ny[list[i][j]];
//      kny[list[i][j]]-= djU*dy + dcU*ny[i];

      kth[i]         -= diU*(nx[i]*dy-ny[i]*dx)+dcU*(nx[i]*ny[list[i][j]]-ny[i]*nx[list[i][j]]);  
	    kth[list[i][j]]-= djU*(nx[list[i][j]]*dy-ny[list[i][j]]*dx)+dcU*(nx[list[i][j]]*ny[i]-ny[list[i][j]]*nx[i]);

	   U1=4.*pow(e1,eta)*pow(e2,mu)*A;
	  }
        }
    }
 
 return 0; 
}



int eq_motion(double* x, double* y,double* theta, double* vx, double* vy,double* nx,double* ny,double* nx_pr,double* ny_pr, double *omega, double dt, double* kx, double* ky, double*knx, double* kny, double* kth, int Np, double* avK, double Th) {
  double zeta,zeta_r;
  double n[Np];
  double n_pr2[Np];
  double R[Np];
  int k;
  for (k = 0; k < Np; k++) {
    n[k] = 0.0;
    R[k] = 0.0;
    n_pr2[k] = 0.0;
  }
  zeta = 1.;
  zeta_r =  100.;
  *avK=0.0;
  for (k = 0; k < Np; k++) {
    vx[k] += -vx[k] * zeta * dt + kx[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vy[k] += -vy[k] * zeta * dt + ky[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    R[k] = sqrt(2. * zeta_r * Th * dt) * gaussian_rand();
    n_pr2[k] = nx_pr[k]*nx_pr[k] + ny_pr[k]*ny_pr[k];
    nx_pr[k] += -n_pr2[k]*nx[k]*dt - zeta_r*nx_pr[k]*dt - kth[k]*ny[k]*dt - ny[k]*R[k];//-nx_pr[k] * zeta_r * dt + knx[k] * dt + ny[k] * R[k];
    ny_pr[k] += -n_pr2[k]*ny[k]*dt - zeta_r*ny_pr[k]*dt + kth[k]*nx[k]*dt + nx[k]*R[k];//-ny_pr[k] * zeta_r * dt + kny[k] * dt - nx[k] * R[k];

    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;

    nx[k] += nx_pr[k] * dt;
    ny[k] += ny_pr[k] * dt;

   n[k] = sqrt(nx[k]*nx[k] + ny[k]*ny[k]);

   nx[k] = nx[k]/ n[k];
   ny[k] = ny[k]/ n[k];

  // cout << nx_pr[1] << "    " << ny_pr[1] <<  "   " << n_pr2[1] << " " << (nx_pr[k] * nx_pr[k] + ny_pr[k] * ny_pr[k])/(n[k]*n[k]) << endl;

    *avK += vx[k]*vx[k] + vy[k]*vy[k] + (nx_pr[k]*nx_pr[k] + ny_pr[k]*ny_pr[k])/(n[k]*n[k]);

    // cout << nx_pr[k] * nx_pr[k] + ny_pr[k] * ny_pr[k] << endl;
  }
  *avK = *avK / Np / 2.0/ 1.5;
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
  int Np = 10000;int count_th=200;
  //double r1=1.0, r2=1.0;
  double disp_th = sqrt(320);
  double dt =0.001;//  //parameters;
  double Th;
  double temp = 1.0;
  double phi = 0.3;
  double dphi = 0.001;

  double kappa=3.; // shape anisotropy parameter
  double kappa_pr=3.; // energy anisotropy parameter
  double chi = (kappa*kappa-1.)/(kappa*kappa+1.); //atof(argv[1]);  
 
  double eta = 1.0;
  double mu = 2.0;
  double chi_pr = (pow(kappa_pr,1./mu)-1.)/(pow(kappa_pr,1./mu)+1.) ; 

  double timer;
  double chi0=0.2;

  double RCHK=4.0; // cut off +1
  double L;
  int M;
   L = (1.0/2.0)*sqrt(pow(2.,1./3.)*pi*kappa*Np/phi);
   M=(int)(L/RCHK);

  //double L = (r1+r2)/4.0*sqrt(pow(2.,1./3.)*pi*kappa*Np/phi);
  //int    M=(int)(L/RCHK);
  //cout << "L=" << L <<" "<< "M="<<M <<endl;
  
  double* x, * y,* x0, * y0, * vx, * vy, * a, * kx, * ky,* knx,* kny,*kth,* nx,* ny,* nx_pr, * ny_pr,*x_update,*y_update,*theta,*omega;
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
  knx = new double[Np];
  kny = new double[Np];
  kth = new double[Np];
  nx = new double[Np];
  ny = new double[Np];
  nx_pr = new double[Np];
  ny_pr = new double[Np];
  theta = new double[Np];
  omega = new double[Np];

  char filename[128];
  avU0=0.0;

  //cout << "L=" << L <<" "<< "M="<<M << " "<< "phi = " << phi <<endl;
  // ini_coord_square(x,y,theta,a,L,Np);

  ini_hex(x,y,nx,ny,theta,a,L,Np);
  
  ini(vx, vy, nx_pr,ny_pr, Np,omega);  

  
  sprintf(filename,"%s/energy_time.txt",argv[1]);
  ofstream file;
  file.open(filename);
  
  update(L,Np,x,y,M,RCHK,list);
  disp_ave = 0.0;

  count=0;
  timer=0.0;
  
  copy(x_update,y_update,x,y,Np,x_corr,y_corr);
  copy(x0,y0,x,y,Np,x_corr,y_corr);

   for(;;){ // infinite loop

   for(int count=0;count<1000;count++){
    calc_force(x, y, nx,ny, L, Np, a, kx, ky,knx,kny,kth,list,theta,chi,chi_pr,mu,eta);    
    eq_motion(x, y, theta, vx, vy, nx, ny, nx_pr, ny_pr, omega, dt, kx, ky,knx, kny,kth, Np, &avK0,temp);
    com_correction(x,y,&x_corr,&y_corr,Np, L);
    p_bound(x, y, Np, L);
    cout <<"disp="<< disp_ave <<"x="<<x[1]<<" "<< "temp=" << avK0 << "L=" << L <<" "<< "M="<<M << " "<< "phi = " << phi << "count = " << count <<endl;
    count++;
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
    
    affine_coord(x, y, phi, Np,dphi);  
       phi += dphi;
     //phi += 0.001*exp(-phi);   
      L = (1.0/2.0)*sqrt(pow(2.,1./3.)*pi*kappa*Np/phi);
      M=(int)(L/RCHK); 
    
    if(phi>=1.0)
    break;

    }


    disp_ave = 0.0;
    int count2 = 0;
    for(;;){
    calc_force(x, y,nx, ny, L, Np, a, kx, ky,knx,kny,kth,list,theta,chi,chi_pr,mu,eta);  
    eq_motion(x, y,theta,vx, vy, nx, ny, nx_pr, ny_pr,omega, dt, kx, ky, knx, kny,kth, Np, &avK0,temp);
    com_correction(x,y,&x_corr,&y_corr, Np, L);
    p_bound(x, y, Np, L);
    cout <<"disp="<< disp_ave <<"x="<<x[1]<<" "<< "temp=" << avK0 << "L=" << L <<" "<< "M="<<M << " "<< "phi = " << phi << " "<<"count2=" << count2<<endl;
     if(count2 > 100000){
    // if(disp_ave>disp_th2){ 
      output(x,y,nx,ny,x_corr,y_corr,L,theta,a,Np); 
      break;
      }
   ///////auto update////////////////////
    calc_disp(&disp_max,&disp_ave,x,y,x_update,y_update,x0,y0,Np,L,x_corr,y_corr);
     count2++;
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
   
  delete[] knx;
  delete[] kny;
  delete[] nx;
  delete[] ny;
  delete[] nx_pr;
  delete[] ny_pr;
  
  delete[] list;
  delete[] theta;
  delete[] omega;
  return 0;
}