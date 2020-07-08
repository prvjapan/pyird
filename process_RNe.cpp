/* developed by T.Hirano, Jan. 2019 */
/* wrapped using boost by H.Kawahara July 2020 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "nrutil_d.c"
#include "nrutil_d.h"
#include "nr_d.h"
#include "boost/python/numpy.hpp"
namespace py = boost::python;
namespace np = boost::python::numpy;

/***************************************/
void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
	  int ma, double **covar, double *chisq, //changed
	  void (*funcs)(double, double [], int));
void covsrt(double **covar, int ma, int ia[], int mfit);//changed
void gaussj(double **a, int n, double **b, int m);
void sort(unsigned long n, double arr[]);

void funcs(double x,double afunc[],int ma);
/***************************************/

const int maxdata = 100001;

/* check!! */
int npix = 2048;
int nchan = 32;
int nRN = 64;
int nlight = 250;
unsigned long nsort = 63;
int n_poly = 2;//order of polynomial
int gain = 2.8;
int rej = 5.0; //clipping level for the polynomial fit

/* main */
np::ndarray processRN(np::ndarray aar,np::ndarray bar){
  //void processRN(np::ndarray aar,np::ndarray bar){

  using namespace std;

  int ia[20],itemp,nlow,nupp,niter,ktemp,ptemp,qtemp,nused;
  double **flux,**fmask,**RN,**fout,**A;
  /*  double **RN,**fout,**A; */

  double median,y[maxdata],fsort[maxdata],ferr[maxdata],RNmaster[maxdata],scale[maxdata];
  double chisq,model,sigma;
  double **covar,par[100];
  double Amed[maxdata],Asigma[maxdata];
  double bias[maxdata],stripe[maxdata];
  double yused[maxdata],fused[maxdata],ferrused[maxdata];

  unsigned long ntemp;
  int ii=0;

  flux = matrix(1,npix+1,1,npix+1);
  fmask = matrix(1,npix+1,1,npix+1);
  RN = matrix(1,nchan+1,1,nRN+1);
  A = matrix(1,nchan+1,1,nlight+1);
  fout = matrix(1,npix+1,1,npix+1);

  covar = matrix(1,n_poly+1,1,n_poly+1);
  
  cin.precision(13);
  cerr.precision(10);
  cout.precision(10);


  /* read the data */
  double *fluxp = reinterpret_cast<double *>(aar.get_data());
  
  for(int i=1;i<=npix;i++){ //y direction
    for(int j=1;j<=npix;j++){ //x direction
      flux[j][i] = fluxp[ii];
      ii=ii+1;
    }
  }

  /* read the mask */
  double *fmaskp = reinterpret_cast<double *>(bar.get_data());
  ii=0;
  for(int i=1;i<=npix;i++){
    for(int j=1;j<=npix;j++){
      fmask[j][i] = fmaskp[ii];
      ii=ii+1;
    }
  }

  
  for(int i=1;i<20;i++){
    ia[i] = 1;
  }

  
  
  /* compute the A factor */

  for(int j=1;j<=nlight;j++){

    for(int k=5;k<=npix;k++){
      y[k-4] = (double)k;
      fsort[k-4] = flux[j][k];
      //ferr[k-4] = sqrt(fabs(flux[j][k])*gain)/gain;
      ferr[k-4] = 10.0;
    }

    lfit(y,fsort,ferr,npix-4,par,ia,n_poly+1,covar,&chisq,funcs);

    /*
    for(int i=1;i<=npix;i++){ //check the starting pix

      if(nsort<=i && i<=npix-nsort/2){
	for(int k=1;k<=nsort;k++){
	  fsort[k] = flux[j][i-nsort/2-1+k];
	}
      }else if(i<nsort){
	for(int k=1;k<=nsort;k++){
	  fsort[k] = flux[j][i+k-1];
	}
      }else if(npix-nsort/2<i){
	for(int k=1;k<=nsort;k++){
	  fsort[k] = flux[j][npix-k];
	}
      }
      sort(nsort,fsort);
      median = fsort[nsort/2+1];

      
      model = par[1];
      for(int k=1;k<=n_poly;k++){
	model += par[k+1]*pow((double)i,(double)k);
      }

      if(j==100){
	//cout << i << " " << flux[j][i] << " " << median << endl;
	cout << i << " " << flux[j][i] << " " << model << endl;
      }

    }
    */

    sigma = 0.0;
    for(int i=5;i<=npix;i++){
      model = par[1];
      for(int k=1;k<=n_poly;k++){
	model += par[k+1]*pow(y[i-4],(double)k);
      }
      sigma += (fsort[i-4]-model)*(fsort[i-4]-model);
    }
    sigma = sqrt(sigma/((double)npix-4.0));

    
    nused = 0;
    for(int i=5;i<=npix;i++){
      model = par[1];
      for(int k=1;k<=n_poly;k++){
	model += par[k+1]*pow(y[i-4],(double)k);
      }
      if(fabs(fsort[i-4]-model)<rej*sigma){
	nused++;
	yused[nused] = (double)i;
	fused[nused] = fsort[i-4];
	ferrused[nused] = ferr[i-4];
      }
    }

    //cerr << "nused = " << nused << " (npix-4 = 2044)" << endl;
    
    lfit(yused,fused,ferrused,nused,par,ia,n_poly+1,covar,&chisq,funcs);
    
    
    for(int p=1;p<=nchan;p++){
      for(int q=1;q<=nRN;q++){
	itemp = nRN*(p-1) + q; 
	
	model = par[1];
	for(int k=1;k<=n_poly;k++){
	  model += par[k+1]*pow((double)itemp,(double)k);
	}

	RN[p][q] = flux[j][itemp]-model;

      }
    }

    ntemp = nchan/2;
    for(int q=1;q<=nRN;q++){
      for(int p=1;p<=ntemp;p++){
	fsort[p] = RN[p*2][q];
      }
      sort(ntemp,fsort);

      RNmaster[q] = (fsort[ntemp/2]+fsort[ntemp/2+1])*0.5;

      if(j==100){
	//cout << q << " " << RNmaster[q]<< endl;
      }
      
    }

     
    /* compute coefficienta A */
    
    for(int p=1;p<=nchan/2;p++){
      for(int q=1;q<=nRN;q++){
	//itemp = nRN*(2*p-1) + q;
	scale[q] = RN[2*p][q]/RNmaster[q];
      }
      sort((unsigned long)nRN,scale);

      A[2*p][j] = (scale[nRN/2]+scale[nRN/2+1])*0.5;

    }

    for(int p=1;p<=nchan/2;p++){
      for(int q=1;q<=nRN;q++){
	scale[q] = RN[2*p-1][q]/RNmaster[nRN+1-q];
      }
      sort((unsigned long)nRN,scale);

      A[2*p-1][j] = (scale[nRN/2]+scale[nRN/2+1])*0.5;

    }
    
    if(j==100){
      for(int p=1;p<=nchan;p++){
	//cout << p << " " << A[p][j] << endl;
      }
    }
    
  }

  
  for(int p=1;p<=nchan;p++){
    for(int j=1;j<=nlight;j++){
      fsort[j] = A[p][j];
    }
    sort((unsigned long)nlight,fsort);

    nupp = nlight*0.8413;
    nlow = nlight*0.1587;
    Amed[p] = (fsort[nlight/2]+fsort[nlight/2+1])*0.5;
    Asigma[p] = (fsort[nupp]-fsort[nlow])*0.5;
  
    //cout << p << " " << Amed[p] << " " << Asigma[p] << endl;

  }

  


  /* subtract the master RN */

  for(int j=1;j<=npix;j++){

    //for(int k=5;k<=npix;k++){
    itemp = 5;
    while(fmask[j][itemp]!=0.0){
      itemp++;
    }

    /*
    if(j==1400){
      for(int i=1;i<=npix;i++){
	cout << i << " " << flux[j][i] << " " << fmask[j][i] << endl;
      }
    }
    */
    //cerr << j << " " << itemp << endl;


    for(int p=1;p<=nchan;p++){
      for(int q=1;q<=nRN;q++){
	RN[p][q] = 0.0;
      }
    }
    
    niter = 1;
    while(itemp<npix){
      ktemp = 1;
      while(fmask[j][itemp]==0.0 && itemp<npix){
	y[ktemp] = (double)itemp;
	fsort[ktemp] = flux[j][itemp];
	//ferr[ktemp] = sqrt(fabs(flux[j][itemp])*gain)/gain;
	ferr[ktemp] = 10.0;
	ktemp++;
	itemp++;
      }
      ktemp--;
      //cerr << j<< " " << ktemp << " " << itemp << endl;


      if(10<ktemp){
	      
	lfit(y,fsort,ferr,ktemp,par,ia,n_poly+1,covar,&chisq,funcs);

	/*
	if(j==94){
	  for(int i=1;i<=ktemp;i++){
	    model = par[1];
	    for(int k=1;k<=n_poly;k++){
	      model += par[k+1]*pow((double)y[i],(double)k);
	    }
	    cout << y[i] << " " << fsort[i] << " " << model << endl;
	  }
	}
	*/

	sigma = 0.0;
	for(int i=1;i<=ktemp;i++){
	  model = par[1];
	  for(int k=1;k<=n_poly;k++){
	    model += par[k+1]*pow(y[i],(double)k);
	  }
	  sigma += (fsort[i]-model)*(fsort[i]-model);
	}
	sigma = sqrt(sigma/(double)ktemp);

    
	nused = 0;
	for(int i=1;i<=ktemp;i++){
	  model = par[1];
	  for(int k=1;k<=n_poly;k++){
	    model += par[k+1]*pow(y[i],(double)k);
	  }
	  if(fabs(fsort[i]-model)<rej*sigma){
	    nused++;
	    yused[nused] = y[i];
	    fused[nused] = fsort[i];
	    ferrused[nused] = ferr[i];
	  }
	}
    
	lfit(yused,fused,ferrused,nused,par,ia,n_poly+1,covar,&chisq,funcs);

	/*
	if(j==94){
	  for(int i=1;i<=ktemp;i++){
	    model = par[1];
	    for(int k=1;k<=n_poly;k++){
	      model += par[k+1]*pow((double)y[i],(double)k);
	    }
	    cout << y[i] << " " << fsort[i] << " " << model << endl;
	  }
	}
	*/
	
	for(int i=1;i<=ktemp;i++){
	  model = par[1];
	  for(int k=1;k<=n_poly;k++){
	    model += par[k+1]*pow((double)y[i],(double)k);
	  }

	  qtemp = (int)y[i] % nRN;
	  if(qtemp==0) qtemp = nRN;
	  
	  ptemp = ((int)y[i]-qtemp)/nRN + 1;
	  
	  RN[ptemp][qtemp] = fsort[i]-model;


	  if(j==1400 && ptemp%2==1){
	    //cout << ptemp << " " << qtemp << " " << RN[ptemp][qtemp] << endl;
	    //cout << ptemp << " " << 65-qtemp << " " << RN[ptemp][qtemp] << endl;
	  }

	  
	}
	
      }

      
      while(fmask[j][itemp]!=0.0 && itemp<npix){
	itemp++;
      }

      niter++;

    }


    for(int q=1;q<=nRN;q++){

      ktemp = 0;
      for(int p=1;p<=nchan;p++){
	if(RN[p][q]!=0.0 && p%2==0){
	  ktemp++;
	  fsort[ktemp] = RN[p][q]/Amed[p];
	}else if(RN[p][nRN+1-q]!=0.0 && p%2==1){
	  ktemp++;
	  fsort[ktemp] = RN[p][nRN+1-q]/Amed[p];
	}
      }
      sort((unsigned long)ktemp,fsort);

      if(ktemp%2==0){
	RNmaster[q] = (fsort[ktemp/2]+fsort[ktemp/2+1])*0.5;
      }else{
	RNmaster[q] = fsort[ktemp/2+1];
      }


      if(j==1400){
	//cout << q << " " << RNmaster[q] << endl;
      }
      
    }

    for(int p=1;p<=nchan;p++){
      for(int q=1;q<=nRN;q++){
	itemp = nRN*(p-1) + q;
	if(p%2==0){
	  fout[j][itemp] = RNmaster[q]*Amed[p];
	}else{
	  fout[j][itemp] = RNmaster[nRN+1-q]*Amed[p];
	}
      }
    }


    /*
    if(j==2047){
      for(int i=1;i<=npix;i++){
	cout << i << " " << flux[j][i] << " " << fmask[j][i] << " " << fout[j][i] << endl;
      }
    }
    */

    
  }
      

  /* recompute the bias value for each channel */

  for(int p=1;p<=nchan;p++){
    ktemp = 0;
    for(int q=1;q<=nRN;q++){
      itemp = nRN*(p-1) + q;
      for(int j=1;j<=4;j++){ //check!!
	ktemp++;
	fsort[ktemp] =  flux[j][itemp]-fout[j][itemp];
      }
      for(int j=npix-3;j<=npix;j++){
	ktemp++;
	fsort[ktemp] =  flux[j][itemp]-fout[j][itemp];
      }
    }
    sort((unsigned long)ktemp,fsort);

    bias[p] = (fsort[ktemp/2]+fsort[ktemp/2+1])*0.5;

    nupp = ktemp*0.8413;
    nlow = ktemp*0.1587;
    //cout << p << " " << bias[p] << " " << (fsort[nupp]-fsort[nlow])*0.5 <<  endl;
    
  }


  /* compute the vertical stripe for each x */

  for(int j=1;j<=npix;j++){
    ktemp = 0;
    for(int i=1;i<=4;i++){ 
      ktemp++;
      fsort[ktemp] = flux[j][i]-fout[j][i]-bias[1];
      //if(j==150) cerr << fsort[ktemp] << endl;
    }
    for(int i=npix-3;i<=npix;i++){ //check!!
      ktemp++;
      fsort[ktemp] = flux[j][i]-fout[j][i]-bias[1];
      //if(j==150) cerr << fsort[ktemp] << endl;
    }
    
    sort((unsigned long)ktemp,fsort);
      
    stripe[j] = (fsort[ktemp/2]+fsort[ktemp/2+1])*0.5;

    //cout << j << " " << stripe[j] << endl;
    
  }

    
  
  
  /* output the final image */
  py::tuple shapea = py::make_tuple(npix,npix);
  np::dtype dtype = np::dtype::get_builtin<double>();
  np::ndarray aout = np::zeros(shapea, dtype);
  
  ktemp = 0;
  ii=0;
  for(int i=1;i<=npix;i++){
    qtemp = i % nRN;
    if(qtemp==0) qtemp = nRN;	  
    ptemp = (i-qtemp)/nRN + 1;
    
    for(int j=1;j<=npix;j++){
      //cout << "\t" << flux[j][i]-fout[j][i]-bias[ptemp]-stripe[j];
      aout[j-1][i-1]=flux[j][i]-fout[j][i]-bias[ptemp]-stripe[j];
      ktemp++;
      ii=ii+1;
      //if(ktemp%3==0) cout << endl;
    }
    
  }

  return aout;
  

}


/***************************************/
void funcs(double x,double afunc[],int ma)
{
  
  afunc[1] = 1.0;
  
  for(int i=2;i<=ma;i++){
    //afunc[i] = pow(x,i-1);
    afunc[i] = afunc[i-1]*x;
  }

  //return 0;

}

/***************************************/
#define NRANSI
void lfit(double x[], double y[], double sig[], int ndat, double a[], 
	  int ia[], int ma, double **covar, double *chisq, //changed
	  void (*funcs)(double, double [], int))
{
  
  void covsrt(double **covar, int ma, int ia[], int mfit);//changed
  void gaussj(double **a, int n, double **b, int m);
  int i,j,k,l,m,mfit=0;
  double ym,wt,sum,sig2i,**beta,*afunc;

  beta=matrix(1,ma,1,1);
  afunc=vector(1,ma);

  for (j=1;j<=ma;j++)
    if (ia[j]) mfit++;
  if (mfit == 0) nrerror("lfit: no parameters to be fitted");
  for (j=1;j<=mfit;j++) {
    for (k=1;k<=mfit;k++) covar[j][k]=0.0;    
    //cerr << "Fuccck!!" << endl;
    beta[j][1]=0.0;
  }

  //cerr << "Fuck!!!" << endl;

  for (i=1;i<=ndat;i++) {
    (*funcs)(x[i],afunc,ma);
    ym=y[i];
    if (mfit < ma) {
      for (j=1;j<=ma;j++)
	if (!ia[j]) ym -= a[j]*afunc[j];
    }
    sig2i=1.0/SQR(sig[i]);
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
	wt=afunc[l]*sig2i;
	for (j++,k=0,m=1;m<=l;m++)
	  if (ia[m]) covar[j][++k] += wt*afunc[m];
	beta[j][1] += ym*wt;
      }
    }
  }
  for (j=2;j<=mfit;j++)
    for (k=1;k<j;k++)
      covar[k][j]=covar[j][k];
  gaussj(covar,mfit,beta,1);
  for (j=0,l=1;l<=ma;l++)
    if (ia[l]) a[l]=beta[++j][1];
  *chisq=0.0;
  for (i=1;i<=ndat;i++) {
    (*funcs)(x[i],afunc,ma);
    for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
    *chisq += SQR((y[i]-sum)/sig[i]);
  }
  covsrt(covar,ma,ia,mfit);
  free_vector(afunc,1,ma);
  free_matrix(beta,1,ma,1,1);
}
#undef NRANSI


/***************************************/

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
void covsrt(double **covar, int ma, int ia[], int mfit)//changed
{
  int i,j,k;
  double swap;
  
  for (i=mfit+1;i<=ma;i++)
    for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit;
  for (j=ma;j>=1;j--) {
    if (ia[j]) {
      for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			    for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
						  k--;
    }
  }
}
#undef SWAP



/***************************************/
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NRANSI
void gaussj(double **a, int n, double **b, int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;
  
  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for (j=1;j<=n;j++) ipiv[j]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j] != 1)
	for (k=1;k<=n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			   for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
						}
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI


/***************************************/

#define NRANSI
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort(unsigned long n, double arr[])
{
  unsigned long i,ir=n,j,k,l=1;
  int jstack=0,*istack;
  double a,temp;
  
  istack=ivector(1,NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	for (i=j-1;i>=1;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    }
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	  }
      if (arr[l+1] > arr[l]) {
	SWAP(arr[l+1],arr[l])
	  }
      i=l+1;
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
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
#undef NRANSI


BOOST_PYTHON_MODULE(process_RNe) {
  Py_Initialize();
  np::initialize();
  /*  py::def("processRNt", processRNt); */
  py::def("processRN", processRN);

}
