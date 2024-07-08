#include <QTime>
#include <QtMath>
#include <QString>

#if !defined(constant_h)
#include "L_Init/Constant.h"
#endif

#include "o_error_function.h"



//extern char* QSstrcpy(char *dest , QString *qs);

double Accident(double Average, double sigma);
double AccidentAsym(double Peak, double sigmaL, double sigmaR );
double AccidentCut(double Average, double sigma, double cut);
void InitRandom(void);
double MyRandom(void);

double FindValueUniversal(double x, double *array, int dimension, bool LogFlag=false);  // siggest array[dimension][2]
double ErrorFunction(double x);
double ErrorFunctionC(double x);
double find_parabola(double x1, double y1, double x2, double y2, double x3, double y3, double x);
double find_line(double x1, double y1, double x2, double y2, double x);


char *GetNextSymbol(char *s);
char *GetNextDelimeter(char *s);
char *GetSpecialDelimeter(char *s, char D);
char *GetIntFromString(int &V, char *s);
char *GetShortIntFromString(short int &V, char *s);
char *GetDoubleFromString(double &V, char *s);
char *GetFloatFromString(float &V, char *s);
char *GetSpecialStrFromString(char *ptr, char *base, bool make_zero=true);
char *GetSubstringFromString(char *dest, char *base, int maxsize);
char *GetQStrFromString(QString &V, char *ss);

int  calculate_bucks(char* str);

double RootSumSquared(double x1, double x2=0, double x3=0, double x4=0);
double AddToRSS(double x, double sum);

char* eos(char *b);
char* eos2(char *b);
double round_decimal(double v, int sign);
void swap(double& x1, double& x2);
void swapLise(double& x1, double& x2) {swap(x1,x2);};
void swap(int& x1, int& x2);
void swap(short int& x1, short int& x2);
double pow2(double par);
double pow3(double par);
bool set_limits(double x, double &xmin, double &xmax);
extern double pow_int(double x, int y);

double thetaRandom();              // return a value in the range 0 to pi
double phiRandom();                 // return a value in the range 0 to 2pi
double EllipseRandom(int mode, double R, double &x, double &y, double dR);
double XYZsphereRandom(int mode, double R, double &x, double &y, double &z, double dR=0);
// mode      0 - sphere,  1- ball, 2 - sphere with radius error

void RectangleRandom(int mode, double xmin, double xmax,  double ymin, double ymax,
                     double &x, double &y, double dx=0, double dy=0);

double SumPolynom(int order, double X, double *coef);

int isign(int k);
int MakeIntRange(int V, int minV, int MaxV);
double MakeDoubleRange(double V, double minV, double maxV);


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


//==========================================================================
//=======================================================================================
char* eos(char *b)
{ char *s; s=b;

  while(1)  {
      if(*s==0 || *s=='\n' || *s==';' ) break;
      s++;
    }

  return s; }
//==========================================================================
//=======================================================================================
char* eos2(char *b)
{ char *s; s=b;

  while(1)  {
      if(*s==0 || *s=='\n' || *s==';' || *s==':' ) break;
      s++;
    }

  return s; }
//=================================================================================
double round_decimal(double v, int sign)
{
  int i, div=1;
  double shift=0.5;
  long l;

  for(i=0; i<sign; i++)  { shift/=10; div*=10;}

  l=(v+shift)*div;
  v=double(l)/div;
  return v;
}
//=================================================================================
/*
#include <dos.h>

static union REGS rg;

int get_char()
{
 int c;
 while (1) {
  rg.h.ah = 1;
  int86( 0x16, &rg, &rg);
  if( rg.x.flags & 0x40) {
 int86( 0x28, &rg, &rg);
 continue;
 }
  rg.h.ah = 0;
  int86( 0x16, &rg, &rg);
  if( rg.h.al ==0)
 c = rg.h.ah | 128;
  else
 c = rg.h.al;

  break;
  }
 return c;
}

*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

#define PowerRand 28

/*static unsigned long a_rand=1664525L;
static unsigned long c_rand=1013904223L;
static unsigned long jflone=0x3f800000;
static unsigned long jflmsk=0x007fffff;
*/
/* for VAX
static unsigned long jflone=0x00004080;
static unsigned long jflmsk=0xffff007f;
*/

/*static unsigned long a_rand=9301L;
static unsigned long c_rand=49297L;
static unsigned long m_rand=233280L;

static unsigned long a_rand2=8121L;
static unsigned long c_rand2=28411L;
static unsigned long m_rand2=134456L;

//unsigned long m_limit=2<<PowerRand;
static unsigned long idum2=7L;
*/
#define gdIA   16807
#define gdIM   2147483647
//#define gdAM   (1.0/gdIM)
#define gdIQ   127773
#define gdIR   2836
#define gdMASK 123459876

static long idum =13L;
static long double gdAM=1.0/gdIM;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double Accident(double Average, double sigma)
{
  if(sigma<=0) return Average;         // 04/15/2013

  double Accident=0;
  for(int i=0; i<12; i++) Accident+=MyRandom();

  Accident = Average + sigma*(Accident-6.);

  return  Accident;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double AccidentAsym(double Peak, double sigmaL, double sigmaR)
{
  sigmaL = qFabs(sigmaL); sigmaR = qFabs(sigmaR);
  double SumSig = sigmaL + sigmaR;
  if(SumSig < 1e-90) return Peak;

  bool Left = MyRandom()*SumSig < sigmaL;
  double sigma = Left ? -sigmaL : sigmaR;

  double Accident=0;
  for(int i=0; i<12; i++) Accident+=MyRandom();


  Accident = Peak + sigma*qFabs(Accident-6.);

  return  Accident;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double AccidentCut(double Average, double sigma, double cut)
{
  sigma = qFabs(sigma);
  cut   = qFabs(cut  );

  if(cut <= 0.1*sigma ) cut=sigma;

  cut/=sigma;

START:

  double Accident=0;
  for(int i=0; i<12; i++) Accident+=MyRandom();

  double shift = Accident-6.;
  if(qFabs(shift)>cut) goto START;

  Accident = Average + sigma*shift;

  return  Accident;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void InitRandom(void)
{
  QTime timep = QTime::currentTime();


  idum =  timep.msec() + (60*(60*timep.hour() + timep.minute())+timep.second());

  /*int zz=2<<PowerA;

F1:  a_rand=random(zz);
if(a_rand%8==5 || m_limit%a_rand==0) goto F1;

F2: c_rand=random(zz);
if(c_rand%8==5 || m_limit%c_rand==0) goto F2;

m_last=random(zz)*m_limit/zz+random(zz);
MyRandom();*/
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double MyRandom(void)
{
  //double ran, ran2;
  //unsigned long m_rand=(a_rand*m_last+c_rand)%m_limit;
  //	  	  m_last=m_rand;

  //long double t1=m_rand;
  //long double t2=m_limit;

  //return t1/t2;
  /*int L1=rand()%10000;
int L2=rand()%10000;
double t1=L1;
double t2=double(L2)/10000.;
t1+=t2;
t1/=10000.;
//if(t1 >= 1. || t1 <= 0)
        t1+=0.00000000001;
return t1;*/
  /*unsigned long itemp;

idum = a_rand*idum+c_rand;
itemp = jflone | (jflmsk & idum);
rand = (*(float *)&itemp)-1.0;*/

  /*idum = (a_rand*idum+c_rand)%m_rand;
ran  = (double)idum /(double)m_rand;

idum2 = (a_rand2*idum2+c_rand2)%m_rand2;
ran2  = (double)idum2 /(double)m_rand2;

ran=(ran+ran2)/2.;*/


  /*
Minimal" random number generator of Park and Miller. Returns a uniform random deviate
between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
to initialize the sequence; idum must not be altered between calls for successive deviates in
a sequence.*/


  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  long k;
  static long sw=0;
  double ans;

  if(sw++ > 10000000)      // v.9.8.104     // v.9.8.109  this value has been decreased by ten
    {
      InitRandom();
      sw=0;
    }

  idum ^= gdMASK;  //XORing with MASK allows use of zero and other
  k= idum/gdIQ;    //simple bit patterns for idum.

  idum=gdIA*(idum-k*gdIQ)-gdIR*k; // Compute idum=(IA*idum) % IM without overows by Schrage's method.
  if ( idum < 0)
    idum += gdIM;
  ans=gdAM*idum; 			  // Convert idum to a floating result.
  idum ^= gdMASK; 			  //  Unmask before return.
  return ans;
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  /*#define num 32760
double f=random(num);
double ran = f/double(num);
return ran;
  */
}

//===========================================================
void swap(double& x1, double& x2)
{
  double swap=x1;
  x1=x2;
  x2=swap;
}
//===========================================================
void swap(int& x1, int& x2)
{
  int swap=x1;
  x1=x2;
  x2=swap;
}
//===========================================================
void swap(short int& x1,short  int& x2)
{
  short int swap=x1;
  x1=x2;
  x2=swap;
}
/*//===========================================================
void swap(distribution *d1, distribution *d2)
{
distribution *swap = d1;
d1 = d2;
d2 = swap;
}
//===========================================================
void swap(distribution2 *d1, distribution2 *d2)
{

distribution2 *swap = d1;
d1 = d2;
d2 = swap;
} */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double FindValueUniversal(double x, double array[][2], int dimension, bool LogFlag)  // siggest array[dimension][2]
{
  extern double find_parabola(double x1, double y1, double x2, double y2, double x3, double y3, double x);
  extern double find_line(double x1, double y1, double x2, double y2, double x);


  double r=BadValue;
  int i=-1;

  do i++;
  while (array[i][0] < x && i < dimension-1);

  if(i==dimension-1 || i==0) {
      return r;
    }

  double ym0, ym1, ym2=1;
  ym0 = array[i  ][1];
  ym1 = array[i-1][1];
  if(i>1) ym2 = array[i-2][1];

  if(LogFlag && ym1>0 && ym0>0 && ym2>0)
    {
      ym0 = log(ym0);
      ym1 = log(ym1);
      ym2 = log(ym2);
    }
  else  LogFlag=false;

  if( i > 1 ) r=find_parabola(array[i-2][0],ym2,
      array[i-1][0],ym1,
      array[i  ][0],ym0,
      x);

  else 		r=find_line(    array[i-1][0],ym1,
      array[i  ][0],ym0,
      x);

  if(LogFlag) r =exp(r);

  return r;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double  ErrorFunction(double x)
{
  double value;
  double sign_save=1;

  if(x < 0) { sign_save=-1; x=-x; };

  if(x > max_erf_value) return sign_save;

  if(x < local_erf_value)
    {
      value = FindValueUniversal(x,t,ef_rows,false);
    }
  else  {
      value = FindValueUniversal(x,tc,efc_rows,true);
      value = 1. - value;
    }

  value*=sign_save;
  return value;


  //#define  local_erf_value  4.1
  //#define  max_erf_value    11.9
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double  ErrorFunctionC(double x)
{
  double value;
  double sign_save=1;

  if(x < 0) { sign_save=-1; x=-x; };

  if(x > max_erf_value) return sign_save;

  if(x < local_erf_value)
    {
      value = FindValueUniversal(x,t,ef_rows,false);
      if(sign_save<0) value = -value;
      value = 1. - value;
    }
  else  {
      value = FindValueUniversal(x,tc,efc_rows,true);
      if(sign_save<0) value = 2. - value;
    }


  return value;


  //#define  local_erf_value  4.1
  //#define  max_erf_value    11.9


}
//===========================================================
// ================================================================================
double find_parabola(double x1, double y1, double x2, double y2, double x3, double y3, double x)
{
  if(x==x1) return y1;
  if(x==x2) return y2;
  if(x==x3) return y3;

  if(y1==y2 && y2==y3)
    return y1;

  //---------------------------
  if(x1 > x2) { swap(x1,x2); swap(y1,y2); }
  if(x2 > x3) { swap(x2,x3); swap(y2,y3); }
  if(x1 > x2) { swap(x1,x2); swap(y1,y2); }

  if(x1 == x2)
    { y2 = (y1+y2)/2.; return find_line(x2, y2, x3, y3, x);}
  if(x3 == x2)
    { y2 = (y3+y2)/2.; return find_line(x2, y2, x1, y1, x);}

  if(y1==y2)
    {
      if( x1 < x && x < x2 ) return y1;
      if( x2 < x && x < x3 ) return find_line(x2, y2, x3, y3, x);
    }
  if(y2==y3)
    {
      if( x2 < x && x < x3 ) return y2;
      if( x1 < x && x < x2 ) return find_line(x2, y2, x1, y1, x);
    }


  double a11=x1*x1;
  double a12=x1;
  double a21=x2*x2;
  double a22=x2;
  double a31=x3*x3;
  double a32=x3;
  double an1=x*x;
  double an2=x;

  double det = a11*(a22-a32) - a21*(a12-a32) + a31*(a12-a22);

  if(det==0) return find_line(x2,y2,x3,y3,x);

  double expr1 = an1 * ( y1 * (a22-a32)  -
                         y2 * (a12-a32)  +
                         y3 * (a12-a22 )   );

  double expr2 =-an2 * ( y1 * (a21-a31)  -
                         y2 * (a11-a31)  +
                         y3 * (a11-a21 )   );

  double expr3 =       ( y1 * (a21*a32-a31*a22)  -
                         y2 * (a11*a32-a31*a12)  +
                         y3 * (a11*a22-a21*a12)   );

  return (expr1 + expr2 + expr3) / det;
}


// ================================================================================
double find_line(double x1, double y1, double x2, double y2, double x)
{

  double det = x1-x2;

  if(det==0) return y1;

  double expr =  x * (y1-y2) -  y1*x2 + y2*x1;
  return expr / det;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetNextSymbol(char *s)
{
  if(!s) return nullptr;
  char c;

  while(true) {
      c=*s;
      if(c=='\0' || c=='\n' || c=='\r') {return nullptr;}

      if(!( c==' ' || c=='\t' || c==',' || c=='|' || c=='(' ||
            c==')' || c=='{'  || c=='}' || c=='[' || c==']'    )) break;

      s++;
    }
  return s;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetNextDelimeter(char *ss)
{
  if(!ss) return nullptr;
  char c;
  char *s = ss;

  while(true)
    {
      c=*s;
      if(c=='\0' || c=='\n' || c=='\r')   {return nullptr;}
      if( c==' ' || c=='\t' || c==',' || c=='|' || c=='(' ||
          c==')' || c=='{'  || c=='}' || c=='[' || c==']'    ) break;
      s++;
    }

  return s;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetSpecialDelimeter(char *s, char D)
{
  if(!s) return nullptr;
  char c;

  while(true) {
      c=*s;
      if(c=='\0' || c=='\n' || c=='\r')   {return nullptr;}
      if( c==D ) break;
      s++;
    }
  return s;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetIntFromString(int &V, char *ss)
{
  if(!ss) return NULL;

  char *s=GetNextSymbol(ss);

  if(!s) {
      V = BadValue;
      return NULL;
    }

  if(*s > 'A' && *s< 'z') V = BadValue;
  else 			V = atoi(s);
  return GetNextDelimeter(s);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetShortIntFromString(short int &V, char *s)
{
  if(!s) return NULL;
  s=GetNextSymbol(s); if(!s) {V=0;  return NULL; }
  if(*s > 'A' && *s< 'z') V = BadValue;
  else 				V = (short int)atoi(s);
  return GetNextDelimeter(s);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetDoubleFromString(double &V, char *s)
{
  if(!s) return NULL;
  s=GetNextSymbol(s); if(!s){V=BadValue; return NULL;}
  if(*s > 'A' && *s< 'z') V  = BadValue;
  else 				V=atof(s);

  return GetNextDelimeter(s);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetFloatFromString(float &V, char *s)
{
  if(!s) return NULL;
  s=GetNextSymbol(s); if(!s){V=BadValue; return NULL;}
  if(*s > 'A' && *s< 'z') V  = BadValue;
  else 				V=atof(s);

  return GetNextDelimeter(s);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetSpecialStrFromString(char *dest, char *base, bool make_zero)    /// sizeof[dest]=100
{
  if(!base) return NULL;



  char *ptr = strchr(base, '\"');   if(!ptr) return ptr; 	ptr++;

  if(!make_zero) return ptr;


  char *ss = strchr(ptr ,'\"');
  if(!ss) return ss;

  *ss = '\0';  ss++;

  if(strlen(ptr) < 99) strcpy(dest,ptr);
  return ss;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetSubstringFromString(char *dest, char *base, int maxsize)
{
  if(!base) return nullptr;


  char *ptr = base;

  //---------------------------

  while(1)
    {
      char c = *ptr;

      if(c=='\0' || c=='\n' || c=='\r')   {return nullptr;}

      if( c==' ' || c=='\t' || c==',' || c=='|' || c=='(' ||
          c==')' || c=='{'  || c=='}' || c=='[' || c==']'    ) ptr++;
      else break;
    }

  //---------------------------
  int k=0;

  while(1)
    {
      char c = *ptr;

      if( c=='\0'|| c=='\n' || c=='\r' ||
          c==' ' || c=='\t' || c==','  || c=='|' || c=='(' ||
          c==')' || c=='{'  || c=='}'  || c=='[' || c==']'    ) break;

      dest[  k]   = c;
      dest[++k] = '\0';

      if(k-1 >= maxsize) break;

      ptr++;
    }


  return ptr;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *GetQStrFromString(QString &V, char *ss)
{
  V="";

  if(!ss) return nullptr;

  char *ptr = ss;

  //---------------------------

  while(1)
    {
      char c = *ptr;

      if(c=='\0' || c=='\n' || c=='\r')   {return nullptr;}

      if( c==' ' || c=='\t' || c==',' || c=='|' || c=='(' ||
          c==')' || c=='{'  || c=='}' || c=='[' || c==']'    ) ptr++;
      else break;
    }

  //---------------------------

  while(1)
    {
      char c = *ptr;

      if( c=='\0'|| c=='\n' || c=='\r' ||
          c==' ' || c=='\t' || c==','  || c=='|' || c=='(' ||
          c==')' || c=='{'  || c=='}'  || c=='[' || c==']'    ) break;

      V.append(c);

      ptr++;
    }


  return ptr;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int  calculate_bucks(char* str)
{
  int space=0;
  int bucks=0;
  char *s = str;
  while(*s)
    {
      if(*s=='$') bucks++;
      else if(*s==' ') space++;
      s++;
    }

  return bucks+space/3;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
double RootSumSquared2(double x1, double x2)
{
  double sum = 0;

  sum=AddToRSS(x1,sum);
  sum=AddToRSS(x2,sum);

  if (sum <=0 ) return 0;
  return sqrt(sum);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
double RootSumSquared(double x1, double x2, double x3, double x4)
{
  double sum = 0;

  sum=AddToRSS(x1,sum);
  sum=AddToRSS(x2,sum);
  sum=AddToRSS(x3,sum);
  sum=AddToRSS(x4,sum);

  //sum =  x1*x1+x2*x2+x3*x3+x4*x4;   // can't find error !

  if (sum <=0 ) return 0;
  return sqrt(sum);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=

double AddToRSS(double x, double sum)
{
  if(x==0) return sum;
  if(qFabs(x)> 1e150) x=1e150;
  sum += x*x;
  return sum;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow2(double par) {return par*par;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow3(double par) {return par*par*par;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

bool set_limits(double x, double &xmin, double &xmax)
{
  if(x < xmin ) xmin = x;
  if(x > xmax ) xmax = x;
  if(xmin >= xmax ) return false;
  return true;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double thetaRandom()              // return a value in the range 0 to pi
{
  double  v = 2. * MyRandom() - 1.;
  double theta = acos(v);
  return theta;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double phiRandom()              // return a value in the range 0 to 2pi
{
  double phi = 2.*PI * MyRandom();
  return phi;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double XYZsphereRandom(int mode, double R, double &x, double &y, double &z, double dR)
{
  // mode      0 - sphere,  1- ball, 2 - sphere with radius error

  double  r = R;

  if(mode==1) r=R*pow(MyRandom(),1./3.);
  else if(mode==2) r=Accident(R,dR);

  double  theta = thetaRandom();
  double  phi   = phiRandom();

  x = r * cos ( phi ) * sin ( theta );
  y = r * sin ( phi ) * sin ( theta );
  z = r * cos ( theta );

  return r;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double EllipseRandom(int mode, double R, double &x, double &y, double dR)
{
  // mode      0 - hollow ellipse,  1- fill ellipse  2 - ellipse with radius error

  double  r = R;

  if(mode==1) r=R*pow(MyRandom(),0.5);
  else if(mode==2) r=Accident(R,dR);

  double  phi   = phiRandom();

  x = r * cos ( phi );
  y = r * sin ( phi );
  return r;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void RectangleRandom(int mode, double xmin, double xmax,  double ymin, double ymax,
                     double &x, double &y, double dx, double dy)

{
  // mode      0 - hollow rectangle,  1- fill rectangle  2 - rectangle with error

  x=0; y=0;

  double x1  = qMin(xmin,xmax);
  double x2  = qMax(xmin,xmax);
  if(x2==x1) x2 = x1 +0.1;

  double y1  = qMin(ymin,ymax);
  double y2  = qMax(ymin,ymax);
  if(y2==y1) y2 = y1 +0.1;


  double r1 = MyRandom();
  double r2 = MyRandom();


  // 0  x_played   y=y1
  // 1  x_played   y=y2
  // 2  y_played   x=x1
  // 3  y_played   x=x2


  if(mode==0 || mode==2)
    {
      int side = r1*4.;
      switch(side)
        {

        case 0 : x = x1 + (x2-x1)*r2;   y = y1;			break;

        case 1 : x = x1 + (x2-x1)*r2;   y = y2;			break;

        case 2 : x = x1;   		  y = y1 + (y2-y1)*r2; 	break;

        case 3 : x = x2;   		  y = y1 + (y2-y1)*r2; 	break;

        }


      if(mode==2)

        {

          x=Accident(x,dx);

          y=Accident(y,dy);

        }

    }

  else if(mode==1)

    {

      x = x1 + (x2-x1)*r1;

      y = y1 + (y2-y1)*r2;

    }

}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


void		BZERO(double* X, int n){for(int i=0; i< int(n/sizeof(double));i++) *(X+i) = 0.;};

void		BZERON(double* X, int n){for(int i=0;i<n;i++) *(X+i) = 0.;};
void		BZERO(int* X, int n){for(int i=0;i< int(n/sizeof(int));i++) *(X+i) = 0;};
void		BZERON(int* X, int n){for(int i=0;i<n;i++) *(X+i) = 0;};


void		BSET(double* X, double Y, int n){for(int i=0;i < int(n/sizeof(double));i++) *(X+i) = Y;};
void		BSETN(double* X, double Y, int n) {for(int i=0;i<n;i++)*(X+i) = Y;};

void		BSET(int* X, int Y, int n){for(int i=0;i< int(n/sizeof(int));i++) *(X+i) = Y;};
void		BSETN(int* X, int Y, int n){for(int i=0;i<n;i++) *(X+i) = Y;};

void		BCOPY(double* in, double* out, int n){for(int i=0;i<n;i++) *(in+i) = *(out+i);};
void		BCOPY(double* in, const double* out, int n){for(int i=0;i<n;i++) *(in+i) = *(out+i);};
void		BCOPY(int* in, int* out, int n){for(int i=0;i<n;i++) *(in+i) = *(out+i);};

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
double SumPolynom(int order, double X, double *coef)
{
  double somme=0;
  double Xpower=1;

  for(int i=0; i<=order; i++)
    {
      somme += coef[i] * Xpower;
      Xpower *=  X;
    }
  return somme;
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
int isign(int k)
{
  if(k>0)    return 1;
  if(k<0)    return -1;
  return 0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int MakeIntRange(int V, int minV, int maxV)
{
  V = qMax(minV,V);
  V = qMin(maxV,V);
  return V;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double MakeDoubleRange(double V, double minV, double maxV)
{
  V = qMax(minV,V);
  V = qMin(maxV,V);
  return V;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//extern   double find_delta_d3a(double x);
//extern   double find_delta_d3s(double x);

//for (int i=-100; i<100; i++)
//  {
//    qDebug() << i << find_delta_d3a(i) << find_delta_d3s(i) ;
//  }

//typedef std::pair<double,double> dpair;
//struct dpair2
//{
//double x;
//double delta;
//};

//typedef std::vector<dpair2> dvector;


//const dvector db3_p1a_delta=
//		//anode
//	{{-87.83, 1.195},
//	{-76.8, 0.704},
//	{-64.4, 0.489},
//	{-50.4, 0.246},
//	{-27.6, -0.049},
//	{-10.4, -0.247},
//	{-1.2, -0.307},
//	{9.6, -0.194},
//	{19.2, -0.060},
//	{28.8, -0.027},
//	{33.6, -0.110},
//	{38.8, 0.000},
//	{48.7, 0.128},
//	{58, 0.067},
//	{64.68, 0.251},
//	{72.09, 0.223}};

//const dvector db3_p1s_delta =
//	{{-84.8, 0.729},
//	{-71.36, 0.496},
//	{-50.8, 0.440},
//	{-38.4, 0.325},
//	{-32.8, 0.127},
//	{-25.6, 0.003},
//	{-16.8, -0.050},
//	{-8.8, 0.111},
//	{-3.2, 0.314},
//	{3.6, -0.004},
//	{11.2, -0.435},
//	{20.8, -0.302},
//	{31.2, -0.282},
//	{40, -0.235},
//	{48.98, -0.390},
//	{72.09, -0.291}};

////WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//double find_delta_d3(double x, const dvector& dv);


//double find_delta_d3a(double x) {return find_delta_d3(x,db3_p1a_delta);}
//double find_delta_d3s(double x) {return find_delta_d3(x,db3_p1s_delta);}

//double find_delta_d3(double x, const dvector& dv)
//{
//int size = dv.size();

//int i;
//for(i=0; i < size; i++)
//  if(dv[i].x > x) break;

//double v;

//if(i==0)  //extrapolation
//  {
//  v = find_line(dv[0].x, dv[0].delta, dv[1].x, dv[1].delta, x);
//  }
//else if(i==size)
//  {
//  v = find_line(dv[size-2].x, dv[size-2].delta, dv[size-1].x, dv[size-1].delta, x);
//  }
//else    //extrapolation
//  {
//    if(i!=1) i--;

//    v = find_parabola(
//          dv[i-1].x, dv[i-1].delta,
//          dv[i]  .x, dv[i]  .delta,
//          dv[i+1].x, dv[i+1].delta,
//          x
//        );
//  }

//return v;
//}
