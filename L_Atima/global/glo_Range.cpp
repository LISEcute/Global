#include <QtMath>
#include <QtDebug>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "globallib.h"
extern double pow_int(double par, int power);

double RANGE(double A, double Z, double AT, double ZT, double T);
double ENERGY(double A, double Z, double AT, double ZT, double RG);
void FDATA(double Z, double *F, double &COEFFI, double &COEFF);
double FCORR(double Z);

double polynom(double *F, double parameter, int order);  // 0,1,2  quadratic

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double RANGE(double A, double Z, double AT, double ZT, double T) {
  //
  //     CALCULATES RANGES WITH E. HANELT FIT (E. HANELT, PHD THESIS)
  //
  double Y1a, Y1,Y2,Y3,Y4,RG;
  double FCOR;
  double F[11],COEFFI,COEFF;  // increased 1
  //
  // *********************************************************************
  if( T <= 1.0E-6 ) return 0;

  FDATA( ZT, F, COEFFI, COEFF);

  //Y1a = 1.0 +F[1]*Z+F[2]*pow_int(Z,2)+F[3]*pow_int(Z,3)+F[4]*pow_int(Z,4);
  F[0] = 1;
  Y1 =  polynom(F, Z, 4);
  //qDebug() << " Y1a Y1 " << Y1a << Y1;

  Y2 = F[5]+F[6]*Z;
  Y3 = F[7]+F[8]*Z;
  Y4 = F[9]+F[10]*Z;
  //      Y5 = Y3/(2.0*Y4);

  //     MATTER COEFFICIENT
  COEFF = (AT*COEFFI)/pow(ZT,COEFF);

  //     RANGE
  FCOR=FCORR(Z);
  double L10=log10(T);
  double t= Y1*(Y2+Y3*L10+Y4*pow2(L10));
  RG =(A/pow2(Z))*pow(10.,t)*COEFF * FCOR;
  RG = qMax(RG,0.);

  return RG;
}

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double ENERGY(double A, double Z, double AT, double ZT, double RG)
{
  //
  //     CALCULATES ENERGIES FROM RANGES WITH EH FIT (E. HANELT, PHD THESIS
  //
  double Y1a,Y1,Y2,Y3,Y4,Y5;
  double FCOR,T;
  double F[11],COEFFI,COEFF;  // increased 1
  // *********************************************************************
  if( RG <= 1.0E-9 ) return 0;

  FDATA( ZT, F, COEFFI, COEFF);

  F[0]=1;
  Y1 =  polynom(F, Z, 4);
  //Y1a = 1.0 +F[1]*Z+F[2]*pow_int(Z,2)+F[3]*pow_int(Z,3)+F[4]*pow_int(Z,4);

  //qDebug() << " Y1a Y1 " << Y1a << Y1;

  Y2 = F[5]+F[6]*Z;
  Y3 = F[7]+F[8]*Z;
  Y4 = F[9]+F[10]*Z;
  Y5 = Y3/(2.0*Y4);

  //     MATTER COEFFICIENT
  COEFF = (AT*COEFFI)/pow(ZT,COEFF);
  //     ENERGY
  FCOR=FCORR(Z);

  double t=-Y5-sqrt(pow2(Y5)-Y2/Y4+log10(RG/FCOR/COEFF*pow2(Z)/A)/(Y1*Y4));
  T =pow(10.,t);
  //      if( T < 0. ) T = 0.;

  return T;
}
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FDATA(double Z, double *F, double &COEFFI, double &COEFF) {
  //
  //     PARAMETERS FOR THE RANGE-ENERGY FITS
  //
  //
  // *********************************************************************
  if (Z<=5)
    {
      F[1] = -0.128482E-03;
      F[2] = -0.173612E-05;
      F[3] = 0.889892E-07;
      F[4] = -0.705115E-09;
      F[5] = -0.553492E+00;
      F[6] = 0.912049E-02;
      F[7] = 0.268184E+01;
      F[8] = -0.529303E-02;
      F[9] = -0.210108E+00;
      F[10] = 0.774360E-03;
      COEFFI =pow( 4.0,0.98)/9.012;
      COEFF  = 0.98;
    }
  else if( Z>5 && Z<=9 )
    {
      F[1] = 0.667801E-03;
      F[2] = -0.392137E-05;
      F[3] = 0.136917E-06;
      F[4] = -0.972996E-09;
      F[5] = -0.490202E+00;
      F[6] = 0.751599E-02;
      F[7] = 0.261390E+01;
      F[8] = -0.600822E-02;
      F[9] = -0.199549E+00;
      F[10] = 0.731880E-03;
      COEFFI =pow( 6.0,0.98)/12.011;
      COEFF  = 0.98;
    }
  else if( Z>9 && Z<=32)
    {
      F[1] = -0.668659E-04;
      F[2] = -0.185311E-05;
      F[3] = 0.873192E-07;
      F[4] = -0.690141E-09;
      F[5] = -0.530758E+00;
      F[6] = 0.898953E-02;
      F[7] = 0.268916E+01;
      F[8] = -0.533772E-02;
      F[9] = -0.214131E+00;
      F[10] = 0.773008E-03;
      COEFFI =pow( 13.0,0.90)/26.982;
      COEFF  = 0.90;
    }
  else if( Z>32 && Z<=64) {
      F[1] =  1.23639E-03;
      F[2] = -6.13893E-06;
      F[3] = 1.84116E-07;
      F[4] = -1.20551E-09;
      F[5] = -0.263421E+00;
      F[6] = 6.34349E-03;
      F[7] = 2.61081E+00;
      F[8] = -6.38315E-03;
      F[9] = -0.204813E+00;
      F[10] = 6.63267E-04;
      COEFFI =pow( 50.0,0.88)/118.69;
      COEFF  = 0.88;
    }
  else if( Z>64 && Z<=76)
    {
      F[1] = 0.199249E-04;
      F[2] = -0.227944E-05;
      F[3] = 0.105063E-06;
      F[4] = -0.829122E-09;
      F[5] = -0.325062E+00;
      F[6] = 0.975017E-02;
      F[7] = 0.268814E+01;
      F[8] = -0.607419E-02;
      F[9] = -0.218986E+00;
      F[10] = 0.869283E-03;
      COEFFI =pow( 73.0,0.88)/180.95;
      COEFF  = 0.88;
    }
  else {
      F[1] = -0.375861E-03;
      F[2] = -0.373902E-05;
      F[3] = 0.148861E-06;
      F[4] = -0.112159E-08;
      F[5] = -0.166220E+00;
      F[6] = 0.126920E-01;
      F[7] = 0.259061E+01;
      F[8] = -0.725322E-02;
      F[9] = -0.202004E+00;
      F[10] = 0.117942E-02;
      COEFFI =pow( 82.0,0.80)/207.2;
      COEFF  = 0.80;
    }
}
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double FCORR(double Z)
{
  //
  //     CALCULATES A CORRECTION TO THE RANGE FIT OF ECKHARD HANELT
  //     based on measured dE/dx values from thesis Christoph Scheidenberg
  //     Karl-Heinz Schmidt, 19. 12. 1995
  //
  //
  double F[5] ={0.965735686, 9.79114E-3, 3.17099E-3, -6.71227E-4,2.28409E-5};
  double  R = Z*Z / 1000.;
  double poly = polynom(F,R,4);
  double Fcorr = 1./ poly;

//  double Fcorra = 1./( 0.965735686 +
//               9.79114E-3 * R +
//               3.17099E-3 * pow_int(R,2)-
//               6.71227E-4 * pow_int(R,3)+
//               2.28409E-5 * pow_int(R,4));

//  qDebug() << " Fcorra Fcorr " << Fcorra << Fcorr;

  return Fcorr;
}
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//
double polynom(double *F, double parameter, int order)  // 0,1,2  quadratic
{
  if(order<0)  return 0;

  double res=F[0];
  double Xproduct = 1;

  for(int i=1; i<=order; i++)
    {
      Xproduct *= parameter;
      res +=  Xproduct*F[i];
    }

  return res;
}
