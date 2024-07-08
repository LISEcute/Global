#include <QtMath>
#include <stdio.h>
#include <stdlib.h>
#include "globallib.h"
#include <QtDebug>


extern double polynom(double *F, double parameter, int order);  // 0,1,2  quadratic
double fnpol(double X);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void bipoco1s(double X, double CX, double ETH, double THE, double Z, double *BP)
{
//static int counter=0;
//qDebug() << counter++;

  //     binding and polarization correction
  //     1s binding correction for targets

  double G, Ga, FF, PH;
  double F[]= {1,9,31,98,12,25,4.2,0.515};
  // **********************************************************************

  G = polynom(F,X,7) /pow_int(1.+X,9);
//        Ga = ( 1.      +
//              9.  * X +
//              31. * pow_int( X,2 )+
//              98. * pow_int( X,3 )+
//              12. * pow_int( X,4 )+
//              25. * pow_int( X,5 )+
//              4.2 * pow_int( X,6 )+
//              0.515*pow_int( X,7 )    ) /pow_int(1.+X,9);

 // qDebug() << " Ga G " << Ga << G;

        if( CX > .035 ) FF = fnpol(CX);
  else  if( CX <= 0   ) FF = 0;
  else                  FF = -3.*PI/4.*(2.*log(CX)+1.);

  double t=X*X*THE*sqrt(ETH);

  if(fabs(t)<1e-20) PH=0;
  else              PH=FF/t;

  PH = 0.;
  *BP=2.*(G-PH)/Z/THE;
};

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void bipoco2s(double X, double CX, double ETH, double THE, double Z, double *BP)
{
  //     binding and polarization correction
  //     2s binding correction for targets
  //
  double G, Ga, FF, PH;
  double F[]= {1.,9.,31.,49.,162.,63.,18.,1.97};
  // **********************************************************************

  G = polynom(F,X,7) /pow_int(1.+X,9);

  //
  // **********************************************************************
//  Ga = ( 1. +
//        9. * X +
//        31. *pow_int( X,2 )+
//        49. *pow_int( X,3 )+
//        162.*pow_int( X,4 )+
//        63. *pow_int( X,5 )+
//        18. *pow_int( X,6 )+
//        1.97*pow_int( X,7 )     ) /pow_int( (1. + X),9);

  // qDebug() << " Ga G " << Ga << G;


        if( CX > .035 ) FF = fnpol(CX);
  else  if( CX <= 0   ) FF = 0;
  else                  FF = -3. * PI / 4. * (2.*log(CX)+1.);

  double t=X*X*THE*sqrt(ETH);

  if(fabs(t)<1e-20) PH=0;
  else              PH=FF/t;

  PH = 0.;
  *BP=2.*(G-PH)/Z/THE;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void bipoco2p(double X, double CX, double ETH, double THE, double Z, double *BP)
{
  //
  //     binding and polarization correction
  //     2P binding correction for targets
  //
  double G, Ga, FF, PH;
  //
  // **********************************************************************
  double F[]= {1,10,45,102,331,67,58,7.8,0.888};
  // **********************************************************************

  G = polynom(F,X,8) /pow_int(1.+X,9);

//  Ga = ( 1. +
//        10.   * X +
//        45.   *pow_int( X,2 )+
//        102.  *pow_int( X,3 )+
//        331.  *pow_int( X,4 )+
//        67.   *pow_int( X,5 )+
//        58.   *pow_int( X,6 )+
//        7.8   *pow_int( X,7 )+
//        0.888 *pow_int( X,8 ))     /pow_int( 1.+X,9);

   //qDebug() << " Ga G " << Ga << G;

  if( CX > .035 ) FF = fnpol(CX);
  else  if( CX <= 0   ) FF = 0;
  else                  FF=-3. * PI / 4. * (2.*log(CX)+1.);

  double t=X*X*THE*sqrt(ETH);

  if(fabs(t)<1e-20) PH=0;
  else              PH=FF/t;

  PH = 0.;
  *BP=2.*(G-PH)/Z/THE;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double fnpol(double X)
{
  //     needed for POLARIZATION correction
  return      exp(-2.*X)/(.031+0.21*sqrt(X)+0.005*X-0.069*sqrt(X)+0.324*X*X);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW




