#include <QtMath>
#include <stdio.h>
#include <stdlib.h>
#include "globallib.h"

extern double BSA, BETA;

double ETATH[41], THETA[5];
double          SKIB, BP1S,  P1S,  SKIS, BP1SS, P1SS,
                SU2S, BP2SU, P2SU, SU2P, BP2PU, P2PU,
                SS2S, BP2SS, P2SS, SS2P, BP2PS, P2PS,
                SMIU, SMIS,  KIB,  KCB,  LIS,  MIS, KIS;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ionicro(double ZF, double U1S, double U1SS, double U2S,
             double U2P, double U3S, double U3D, int *IRC2)
{
//
//     subroutine for ionization  cross sections / (ZT^2+ZT)
//
      int I,IRC;
      double  ETH1S, ETH1SS, ETH2P, ETH2S;
      double ETH3D, ETH3S, PP, REL, RELK, RELKS, RTR, SI, SS3D, SS3P;
      double SS3S, SU3D, SU3P, SU3S, TH1S, TH1SS, TH2P, TH2S, TH3D, TH3S;
      double TR1S, TR1SS, TR2PS, TR2PU, TR2SS, TR2SU, TR3SS, TR3SU;
      double CX, X, Z1S, Z2P, Z2S, Z3D, Z3S;
//      common /ETATH/  ETATH, THETA;
// ****************************************************************;
//      double         PI, U0, MC2, S0, ALPH, GA, BSQ, BETA, BSA, D1, D2, GAG, ETA;

// ****************************************************************;
// ****************************************************************;
/*	common /BINPOL/ SKIB, BP1S,  P1S,  SKIS, BP1SS, P1SS,
                        SU2S, BP2SU, P2SU, SU2P, BP2PU, P2PU,
                        SS2S, BP2SS, P2SS, SS2P, BP2PS, P2PS,
                        SMIU, SMIS,  KIB,  KCB,  LIS,  MIS;*/
// ****************************************************************;
//  ko vsem massivam dobavlen 0 v pervuju pozitsiju i razmernost' uvelichena na 1
//
//     THETA K in F tables of Benka and Kropf, ADND tables 22 (1978) 219
      double THETAK[5]={0,0.5334,0.7113,0.9486,1.2649};
//     THETA L in F tables of Benka and Kropf, ADND tables 22 (1978) 219
      double THETAL[5]={0,0.4743,.6325,.8434,1.1247};
//     THETA M in F tables of Johnson et al., ADND tables 24 (1979) 1
      double THETAM[5]={0,0.3,0.4,0.5,0.6};
//     ETA/THSQ K and L in F tables of Benka and Kropf
      double ETATHKL[11]={0,.01,.01334,.01778,.02371,.03162,
                          .04217,.05623,.06494,.07499,.0866};
//     ETA/THSQ M in F tables of Johnson et al.
      double ETATHM[41]={0,
             0.04, 0.045, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.3,
              0.4,  0.5,  0.6,  0.7,   0.8,  0.9, 1.,   1.3, 1.5, 1.7,
              2.0,  2.5,  3.0,  3.5,   4.0, 4.5,  5.,   5.5, 6.0, 6.5,
              7.0,  8.0, 10.0, 13.0, 15.0,  20.0, 30.0, 45.0,70.0, 100};
// ****************************************************************;
      Z1S=ZF-0.3;
      Z2S=ZF-1.7;
      Z2P=ZF-4.15;
      Z3S=ZF-8.8;
      Z3D=ZF-21.5;

      if( ZF <= 0. || U1S <= 0. ) {
                TH1S = 0.;
                ETH1S = 0.;
              }
      else    {
                TH1S=U1S/U0/ZF/ZF;
                ETH1S=BSA/pow2(ZF*TH1S);
              }


      if( Z1S <= 0. || U1SS <= 0. ) {
                TH1SS = 0.;
                ETH1SS = 0.;
              }
      else    {
                TH1SS=U1SS/U0/Z1S/Z1S;
                ETH1SS=BSA/pow2(Z1S*TH1SS);
              }

      if( Z2S <= 0. || U2S <= 0. ) {
                TH2S = 0.;
                ETH2S = 0.;
              }
      else    {
                TH2S=U2S/U0/Z2S/Z2S*4.;
                ETH2S=BSA/pow2(Z2S*TH2S);
              }

      if( Z2P <= 0. || U2P <= 0. ) {
                TH2P = 0.;
	        ETH2P = 0.;
              }
      else    {
                TH2P=U2P/U0/Z2P/Z2P*4.;
                ETH2P=BSA/pow2(Z2P*TH2P);
              }

      if( Z3S <= 0. || U3S <= 0. ) {
                TH3S = 0.;
                ETH3S = 0.;
              }
      else    {
                TH3S=U3S/U0/Z3S/Z3S*9.;
                ETH3S=BSA/pow2(Z3S*TH3S);
              }

      if( Z3D <= 0. || U3D <= 0. ) {
                TH3D = 0.;
                ETH3D = 0.;
              }
      else    {
                TH3D=U3D/U0/Z3D/Z3D*9.;
                ETH3D=BSA/pow2(Z3D*TH3D);
              }
//
// *** Compute ioniz. cross section
//L10:   //continue;


      for(I=1; I<=4; I++)  THETA[I] = THETAK[I];

      int IN ;
      for(IN=1; IN<=10; IN++) {
        ETATH[   IN]=      ETATHKL[IN];
        ETATH[10+IN]=10.  *ETATH[IN];
        ETATH[20+IN]=100. *ETATH[IN];
        ETATH[30+IN]=1000.*ETATH[IN];
      }
//
//     unscreened 1s cross section
//L100:
      IRC = *IRC2 ; // Fehlerflag

      intpoklm(ETH1S, TH1S, ZF,  U1S, "fk0", &SI, &TR1S, &PP, &IRC) ;

      if (IRC != 0 ) {
                 IRC = -1;
                 SKIB = 0.;
                 P1S = 0.;
                 BP1S = 0.;
               }
      else     {
                                //     per electron and per (ZT^2+ZT),contains transverse
                 KIB=SI/2.;
                                //     sq. root of 1/relativity corr.
                 REL=1.+.0000133*BETA*ZF*ZF;
                 RELK=1./REL/REL;
                                //     X and CX needed for binding-polarization correction
          	 SKIB=KIB*RELK;
                 P1S=PP-1.;
                 X=2.*sqrt(ETH1S);
                 CX=1.5/X;
                 bipoco1s( X, CX, ETH1S,  TH1S,  ZF,  &BP1S);
               }

//L200:
//     screened 1s cross section
      intpoklm(ETH1SS,TH1SS,Z1S, U1SS,"fk0",&SI,&TR1SS,&PP,&IRC) ;

      if (IRC != 0 ) {
                 IRC = -2;
                 SKIS = 0.;
                 P1SS = 0.;
                 BP1SS = 0.;
               }
      else     {
                        //     per electron and per (ZT^2+ZT),contains transverse
                 KIS=SI/2.;
                        //     sq. root of 1/relativity corr.
                 REL=1.+.0000133*BETA*Z1S*Z1S;
                 RELKS=1./REL/REL;
                 SKIS=KIS*RELKS;
                 P1SS=PP-1.;
                        //     X and CX needed for binding-polarization correction
                 X=2.*sqrt(ETH1SS);
                 CX=1.5/X;
                 bipoco1s( X, CX, ETH1SS, TH1SS, Z1S, &BP1SS);
               }
//L300:
// *** Compute unscreened 2s and 2p cross sections from 2s data
	for(I=1; I<=4; I++)
                THETA[I] = THETAL[I];


        intpoklm(ETH2S, TH2S, Z2S, U2S, "fl1", &SU2S, &TR2SU, &PP, &IRC) ;

	if (IRC != 0 ) {
                 IRC = -3;
                 SU2S = 0.;
                 P2SU = 0.;
                 BP2SU = 0.;
               }
        else   {
                 P2SU=PP-1.;
                 X=4.*sqrt(ETH2S);
                 CX=3./X;
                 bipoco2s( X, CX, ETH2S,  TH2S,  Z2S, &BP2SU);
               }

//L400:
      intpoklm(ETH2S, TH2S, Z2S, U2S, "fl2", &SU2P, &TR2PU, &PP, &IRC) ;

      if (IRC != 0 ) {
                IRC = -4;
                SU2P = 0.;
                P2PU = 0.;
        	BP2PU = 0.;
              }
       else   {
                P2PU=PP-1.;
                X=4.*sqrt(ETH2S);
                CX=2.5/X;
                        //     no relativity correction for L shell
                bipoco2p(X, CX, ETH2S,  TH2S,  Z2S, &BP2PU);
              }
//
//L500:
// *** Compute screened 2s and 2p cross sections from 2p data.
	intpoklm(ETH2P, TH2P, Z2P, U2P, "fl1", &SS2S, &TR2SS, &PP, &IRC) ;

	if (IRC != 0 ) {
                 IRC = -5;
                 SS2S = 0.;
                 P2SS = 0.;
                 BP2SS = 0.;
               }
        else   {
                 P2SS=PP-1.;
                 X=4.*sqrt(ETH2P);
                 CX=3./X;
                 bipoco2s( X, CX, ETH2P,  TH2P,  Z2P, &BP2SS);
               }
//L600:
      intpoklm(ETH2P, TH2P, Z2P, U2P, "fl2", &SS2P, &TR2PS, &PP, &IRC) ;
      if (IRC != 0 ) {
                 IRC = -6;
                 SS2P = 0.;
                 P2PS = 0.;
                 BP2PS = 0.;
               }
      else     {
                 P2PS=PP-1.;
                 X=4.*sqrt(ETH2P);
                 CX=2.5/X;
                 bipoco2p( X, CX, ETH2P,  TH2P,  Z2P, &BP2PS);
               }
//
//L700:
//     Compute unscreened 3s,3p,3d cross sections from 3s data

      for(I=1; I<=4; I++)      THETA[I] = THETAM[I];

      for(I=1; I<=40; I++)     ETATH[I] = ETATHM[I];

      intpoklm(ETH3S, TH3S, Z3S, U3S, "fm1", &SU3S, &TR3SU, &PP, &IRC);
      intpoklm(ETH3S, TH3S, Z3S, U3S, "fm2", &SU3P, &RTR,   &PP, &IRC);
      intpoklm(ETH3S, TH3S, Z3S, U3S, "fm3", &SU3D, &RTR,   &PP, &IRC);

      if( IRC != 0 ) IRC = -7;

      SMIU=(SU3S+SU3P+SU3D)/18.;        //     per electron


//
//L800:
// *** Compute screened 3s,3p,3d cross sections from 3d data
      intpoklm(ETH3D, TH3D, Z3D, U3D, "fm1", &SS3S, &TR3SS, &PP, &IRC);
      intpoklm(ETH3D, TH3D, Z3D, U3D, "fm2", &SS3P, &RTR,   &PP, &IRC);
      intpoklm(ETH3D, TH3D, Z3D, U3D, "fm3", &SS3D, &RTR,   &PP, &IRC);

      if( IRC != 0 ) IRC = -8;

      SMIS=(SS3S+SS3P+SS3D)/18.;     //     per electron

      *IRC2 = IRC ;
}
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$;
