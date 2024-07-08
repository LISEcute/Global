// http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf.html

#include <QtMath>
#include <stdlib.h>
#define TINY 1.0e-20  // A small number.
#include <QDebug>
//#define pow(x,y)  exp((y)*log(x))

#ifndef _min_max
    #include <QtMath>
#endif



//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ludcmp (double **a, int n, int *indx, double *d);
void ludcmp0(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
void lubksb0(double **a, int n, int *indx, double *b);
void savgol(double *c, int np, int nl, int nr, int ld, int m);

double* vector(int start, int dimension);
int* ivector(int start, int dimension);
double** matrix(int startx, int dimx, int starty, int dimy);
void free_vector(double *vv, int start, int dimension);
void free_ivector(int *vv, int /*start*/, int /*dimension*/);
void free_matrix(double **aa, int startx, int dimx, int /*starty*/, int /*dimy*/);

double* vector0(int dimension);
int* ivector0(int dimension);
double** matrix0(int dimx,  int dimy);
void free_vector0(double *vv);
void free_ivector0(int *vv);
void free_matrix0(double **aa,  int dimx, int dimy);

int pow_int(int x, int y);
double pow_int(double x, int y);
double bessk0(double x);
double bessk1(double x);
double bessi0(double x);
double bessi1(double x);
void polcoe(double *x, double *y, int n, double *cof);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
void nrerror(const char *error_text)
{
qDebug() << "NR_Math: " << error_text;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
/*Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation eected by the partial
pivoting; d is output as  1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix.*/

void ludcmp(double **a, int n, int *indx, double *d)
{
int i,imax=0,j,k;
double big,dum,sum,temp;
double *vv;                // vv stores the implicit scaling of each row.

vv=vector(1,n);
*d=1.0; 				// No row interchanges yet.
for (i=1;i<=n;i++) { 		//Loop over rows to get the implicit scaling information.
	big=0.0;
	for (j=1;j<=n;j++) if ((temp=qFabs(a[i][j])) > big) big=temp;
	if (big == 0.0) big=1; //nrerror("Singular matrix in routine ludcmp");  //No nonzero largest element.
	vv[i]=1.0/big; 								   //	Save the scaling.
	}

for (j=1;j<=n;j++)   			//This is the loop over columns of Crout's method.
      {
	for (i=1;i<j;i++) { 		// This is equation (2.3.12) except for i = j.
		sum=a[i][j];
		for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
		}

	big=0.0; 				// Initialize for the search for largest pivot element.
	for (i=j;i<=n;i++) { 		// This is i = j of equation (2.3.12) and i = j +1: ::N
		sum=a[i][j];		//		of equation (2.3.13).
		for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
		if ( (dum=vv[i]*qFabs(sum)) >= big) { //Is the figure of merit for the pivot better than the best so far?
			big=dum;
			imax=i;
			}
		}

	if (j != imax) { 				// Do we need to interchange rows?
		for (k=1;k<=n;k++) { 		// Yes, do so...
			dum=a[imax][k];
			a[imax][k]=a[j][k];
			a[j][k]=dum;
			}
		*d = -(*d); 			//...and change the parity of d.
		vv[imax]=vv[j]; 			// Also interchange the scale factor.
		}

      indx[j]=imax;
	if (a[j][j] == 0.0) a[j][j]=TINY;
			// If the pivot element is zero the matrix is singular (at least to the precision of the
			// algorithm). For some applications on singular matrices, it is desirable to substitute
			// TINY for zero.
	if (j != n) { 		// Now, finally, divide by the pivot element.
		dum=1.0/(a[j][j]);
		for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	} 						 //	Go back for the next column in the reduction.
free_vector(vv,1,n);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
void ludcmp0(double **a, int n, int *indx, double *d)
{
int i,imax=0,j,k;
double big,dum,sum,temp;
double *vv;                // vv stores the implicit scaling of each row.

vv=vector0(n);
*d=1.0; 				// No row interchanges yet.

for (i=0;i<n;i++)  		//Loop over rows to get the implicit scaling information.
	{
    big=0.0;
	for (j=0;j<n;j++) if ((temp=qFabs(a[i][j])) > big) big=temp;
	if (big == 0.0) {big=1; nrerror("Singular matrix in routine ludcmp");}  //No nonzero largest element.
	vv[i]=1.0/big; 								   //	Save the scaling.
	}

for (j=0;j<n;j++)   			//This is the loop over columns of Crout's method.
      {
	for (i=0;i<j;i++)  		// This is equation (2.3.12) except for i = j.
		{
            sum=a[i][j];
		for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
		}

	big=0.0; 				// Initialize for the search for largest pivot element.
	for (i=j;i<n;i++)  		// This is i = j of equation (2.3.12) and i = j +1: ::N
		{
            sum=a[i][j];		//		of equation (2.3.13).
		for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
		if ( (dum=vv[i]*qFabs(sum)) >= big) { //Is the figure of merit for the pivot better than the best so far?
			big=dum;
			imax=i;
			}
		}

	if (j != imax) { 				// Do we need to interchange rows?
		for (k=0;k<n;k++) { 		// Yes, do so...
			dum=a[imax][k];
			a[imax][k]=a[j][k];
			a[j][k]=dum;
			}
		*d = -(*d); 			//...and change the parity of d.
		vv[imax]=vv[j]; 			// Also interchange the scale factor.
		}

      indx[j]=imax;
	if (a[j][j] == 0.0) a[j][j]=TINY;
			// If the pivot element is zero the matrix is singular (at least to the precision of the
			// algorithm). For some applications on singular matrices, it is desirable to substitute
			// TINY for zero.
	if (j != n-1) { 		// Now, finally, divide by the pivot element.
		dum=1.0/(a[j][j]);
		for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	} 						 //	Go back for the next column in the reduction.
free_vector0(vv);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
/*
Solves the set of n linear equations A x X = B. Herea[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n, andindx are not modified by this routine
and can be left in place for successive calls with dierent right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is ecient for use
in matrix inversion.*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
void lubksb(double **a, int n, int *indx, double *b)
{
int i,ii=0,ip,j;
double sum;

for (i=1;i<=n;i++) {   		// When ii is set to a positive value, it will become the
	ip=indx[i];		     //index of the first nonvanishing element of b. Wenow
	sum=b[ip];		     //do the forward substitution, equation (2.3.6). The
	b[ip]=b[i];		     //only new wrinkle is to unscramble the permutation as we go.
	if (ii) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
	else if (sum) ii=i; 		// A nonzero element was encountered, so from now on we
						// will have to do the sums in the loop above.
      b[i]=sum;
	}
for (i=n;i>=1;i--) {         // Now we do the backsubstitution, equation (2.3.7).
	sum=b[i];
	for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
	b[i]=sum/a[i][i];     // Store a component of the solution vector X.
	}			    //	 All done!
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
void lubksb0(double **a, int n, int *indx, double *b)
{
int i,ii=0,ip,j;
double sum;

for (i=0;i<n;i++) {   		// When ii is set to a positive value, it will become the
	ip=indx[i];		     //index of the first nonvanishing element of b. Wenow
	sum=b[ip];		     //do the forward substitution, equation (2.3.6). The
	b[ip]=b[i];		     //only new wrinkle is to unscramble the permutation as we go.
	if (ii!=0) for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
	else if (sum != 0.0) ii=i+1; 		// A nonzero element was encountered, so from now on we
						      // will have to do the sums in the loop above.
      b[i]=sum;
	}
for (i=n-1;i>=0;i--) {         // Now we do the backsubstitution, equation (2.3.7).
	sum=b[i];
	for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
	b[i]=sum/a[i][i];     // Store a component of the solution vector X.
	}			    //	 All done!
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
/*Returns in c[1..np], in wrap-around order (N.B.!) consistent with the argument respns in
routine convlv, a set of Savitzky-Golay filter coecients. nl is the number of leftward (past)
data points used, while nr is the number of rightward (future) data points, making the total
number of data points used nl +nr +1. ld is the order of the derivative desired (e.g., ld = 0
for smoothed function). m is the order of the smoothing polynomial, also equal to the highest
conserved moment; usual values are m = 2or m = 4. */

void savgol(double *c, int np, int nl, int nr, int ld, int m)
{
int imj,ipj,j,k,kk,mm,*indx;
double d,fac,sum,**a,*b;

if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m) return; //nrerror("bad args in savgol");

indx=ivector(1,m+1);
a=matrix(1,m+1,1,m+1);
/*double aa[3][3];

for(k=0; k<3; k++)
	for(kk=0; kk<3; kk++) aa[k][kk]=(k+1)*10+kk+1;
*/
b=vector(1,m+1);

for (ipj=0;ipj<=(m << 1);ipj++) { 		// Set up the normal equations of the desired
	sum=(ipj ? 0.0 : 1.0);			// least-squares fit.
	for (k=1;k<=nr;k++)
      			sum += pow_int(k,ipj);
	for (k=1;k<=nl;k++)
      			sum += pow_int(-k,ipj);
	mm=qMin(ipj,2*m-ipj);
	for (imj = -mm;imj<=mm;imj+=2)
      			a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
	}

for (j=1;j<=m+1;j++) b[j]=0.0;
b[ld+1]=1.0;

ludcmp(a,m+1,indx,&d); 		     // Solve them: LU decomposition.
lubksb(a,m+1,indx,b);		     // Right-hand side vector is unit vector,
                                   // depending on which derivative we want.
			 		     // Get one row of the inverse matrix.
for (kk=1;kk<=np;kk++) c[kk]=0.0;  // Zero the output array (it may be bigger than
					     // number of coeficients).
for (k = -nl;k<=nr;k++) {
	sum=b[1];			     // Each Savitzky-Golay coeficient is the dot
	fac=1.0; 			     //  product of powers of an integer with the inverse matrix row.
	for (mm=1;mm<=m;mm++)
      		sum += b[mm+1]*(fac *= k);
//	kk=((np-k) % np)+1; 		//Store in wrap-around order.     OLEG TARASOV
//	c[kk]=sum;
	c[k+nl]=sum;
	}
free_vector(b,1,m+1);
free_matrix(a,1,m+1,1,m+1);
free_ivector(indx,1,m+1);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double* vector(int start, int dimension)
{
// C++ - start =0, mathematic - start=1

double *v=new double[start+dimension];
return v;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int* ivector(int start, int dimension)
{
// C++ - start =0, mathematic - start=1

int *v=new int[start+dimension];
return v;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double** matrix(int startx, int dimx, int starty, int dimy)

{
// C++ - start =0, mathematic - start=1
double **aa=new double*[startx+dimx];
for(int i=0; i<startx+dimx; i++) aa[i]=new double[starty+dimy];

/*double **aa=new double*[starty+dimy];
for(int i=0; i<starty+dimy; i++) aa[i]=new double[startx+dimx];
*/
return aa;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void free_vector(double *vv, int /*start*/, int /*dimension*/)
{
delete vv;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void free_ivector(int *vv, int /*start*/, int /*dimension*/)
{
delete vv;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void free_matrix(double **aa, int startx, int dimx, int /*starty*/, int /*dimy*/)
{
for(int i=0; i<startx+dimx; i++) delete aa[i];
delete aa;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double* vector0(int dimension)
{
double *v=new double[dimension];
return v;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int* ivector0(int dimension)
{
int *v=new int[dimension];
return v;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double** matrix0(int dimx, int dimy)

{
double **aa=new double*[dimx];
for(int i=0; i<dimx; i++) aa[i]=new double[dimy];
return aa;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void free_vector0(double *vv)
{
delete vv;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void free_ivector0(int *vv)
{
delete vv;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void free_matrix0(double **aa, int dimx, int /*dimy*/)
{
for(int i=0; i<dimx; i++) delete aa[i];
delete aa;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int pow_int(int x, int y)
{
if(y < 0) return 0; // don't use

int sum=1;
for(int i=1; i<=y; i++) sum*=x;
return sum;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow_int(double x, int y)
{
if(y < 0) return 0; // don't use

double sum=1;
for(int i=1; i<=y; i++) sum*=x;
return sum;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*
M nL nR Sample Savitzky-Golay Coefficients
2 2 2 . 0.086 0.343 0.486 0.343 . 0.086
2 3 1 . 0.143 0.171 0.343 0.371 0.257
2 4 0 0.086 . 0.143 . 0.086 0.257 0.886
2 5 5 . 0.084 0.021 0.103 0.161 0.196 0.207 0.196 0.161 0.103 0.021 . 0.084
4 4 4 0.035 . 0.128 0.070 0.315 0.417 0.315 0.070 . 0.128 0.035
4 5 5 0.042 . 0.105 . 0.023 0.140 0.280 0.333 0.280 0.140 . 0.023 . 0.105 0.042
*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double  bessk0(double x)
{
double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
return ans;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double bessk1(double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
return ans;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double bessi0(double x)
{
	double ax,ans;
	double y;

	if ((ax=qFabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double bessi1(double x)
{
	double ax,ans;
	double y;

	if ((ax=qFabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} else {
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void polcoe(double *x, double *y, int n, double *cof)
{
    int k,j,i;
    double phi,ff,b,*s;

    s=vector(0,n);

    for (i=0; i<n; i++)
    			s[i]=cof[i]=0.0;

    s[n-1] = -x[0];

    for (i=1; i<n; i++)
    	  {
        for (j=n-i-1; j<n-1; j++)
            s[j] -= x[i]*s[j+1];

        s[n-1] -= x[i];
    	  }

    for (j=0; j<n; j++)
    	  {
        phi=n;

        for(k=n-1; k>0; k--)
	            phi=k*s[k]+x[j]*phi;

        ff=y[j]/phi;
        b=1.0;

        for(k=n-1; k>=0; k--)
        	{
            cof[k] += b*ff;
            b=s[k]+x[j]*b;
        	}
        }

    free_vector(s,0,n);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

