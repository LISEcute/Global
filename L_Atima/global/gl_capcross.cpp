#include <QtMath>
#include <stdio.h>
#include <stdlib.h>
#include "globallib.h"

extern double BSA, GA, BETA, BSQ, D1a, D2, GAG, ETA;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void capcross(double ZP, double ZT, double *SC, int *IRC )
{
//
//     subroutines for capt. cross
//
//     Eikonal approximation for relativistic non-radiative capture
//     based on analytical expression of Eichler. For K+L+M target
//     electrons into bare projectile K,L,M shells.
//     NOTATION: SEIK(I) eikonal cross. from all target shells into
//     bare proj. I shell. SEIK(I,K) eikonal cross. from target I to
//     bare proj. K shell.
//
// **********************************************************************
//  dimension 3 -> changed everywhere for 4
      double Z2PFAC[4][4], OBKFAC[4][4], EIKFAC[4][4], EIK[4][4], MAG[4][4];
      double ORB[4][4], CEIK[4][4], EXPEZ, EXZT, EZP, EZP2, PEZ, PM;
      double ZFAC, ZPS, ZTS, Z1, Z2, Z2Q, EF, EI;
      int IL, IM, JP, JT;
// **********************************************************************
// **********************************************************************

      if( ZT < 10. ) IL=(int)(ZT-2.); 
      else           IL=8;
           if( ZT < 19. ) IM=0;
	   else if( ZT < 28. ) IM=(int)(ZT-10.); 
      else                IM=18;      

      for(JP=1; JP<=3; JP++) {     //L10

        ZPS=ZP/((double)(JP));  

        for(JT=1; JT<=3; JT++) {    //L20
                 if( JT == 1 ) ZTS=ZT-0.3;
	    else if( JT == 2 ) ZTS=(ZT-3.0)/2.;
            else if( JT == 3 ) ZTS=(ZT-19.0)/3.;

	    if( JT == 1 && JP == 1 )  Z2PFAC[JP][JT] = 1.;
            else                      Z2PFAC[JP][JT]=1.16;

            if( ZPS <= ZTS) {
                    Z1=ZPS;
                    Z2=ZTS;
                    Z2Q=ZTS*Z2PFAC[JP][JT];
                  }
            else  {
                    Z1=ZTS;
                    Z2=ZPS;
                    Z2Q=ZPS*Z2PFAC[JP][JT];
                  }

          if(Z1<=0 || Z2<=0) {             //  Oleg Tarasov inserted
                        CEIK[JP][JT]=0;
                        continue;
                      }

          EI = sqrt(1.-ALPH*ALPH*Z2*Z2);
          EF = sqrt(1.-ALPH*ALPH*Z1*Z1);
          PM = ETA*(EF/GA-EI)/ALPH/ALPH;
	  EZP= ETA*Z2Q;
          EZP2=EZP*EZP;
          PEZ =PI*EZP;
          EXPEZ=exp(PEZ);
          EXZT=-2.*EZP*atan(-PM/Z2);
          ZFAC=Z1*Z2/(Z2*Z2+PM*PM);
          OBKFAC[JP][JT]=2.8E+07*128.*PI*ETA*ETA*pow_int(ZFAC,5)/GAG/GA/5.;
          EIKFAC[JP][JT]=2.*PEZ*exp(EXZT)/(EXPEZ-1./EXPEZ);
          EIK[JP][JT]=1.+5.*EZP*PM/Z2/4.+5.*EZP2*PM*PM/Z2/Z2/12.+EZP2/6.;
          MAG[JP][JT]=-D2+5.*D2*D2/16.+5.*D2*GAG*Z2Q/Z2/8.+D2*EZP2/4.+5.*D2*D2*EZP2/48.;

		  ORB[JP][JT]=5.*PI*D1a*ALPH*(Z1+Z2)*(1.-D2/2.)/18.-
					  5.*D1a*ALPH*Z2*EZP*(1.-D2/2.)/8.     -
					  5.*PI*D1a*GAG*ALPH*Z1*Z2Q/18./Z2     +
					  5.*PI*D1a*GAG*GAG*ALPH*Z1*Z2Q*Z2Q/28./Z2/Z2-
					  5.*PI*D1a*GAG*ALPH*(Z1+Z2-D2*Z1)*Z2Q/28./Z2;

          CEIK[JP][JT]=OBKFAC[JP][JT]*EIKFAC[JP][JT]*(EIK[JP][JT]+MAG[JP][JT]+ORB[JP][JT]);
        } //L20

        SC[JP]=JP*JP*(2.*CEIK[JP][1]+IL*CEIK[JP][2]+IM*CEIK[JP][3]);
 
      } //L10

  *IRC = 0;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void reccork(double Z, double U, double *SRKBEST, int *IRC )
{
//
//     REC subroutine for K shell
//
	  double B, NU, E, KH, NS, NB, DB, RECF, /*AUNP,*/ SRK, AUSA;
      double SRKSA;
// **********************************************************************
//     start of K REC calculation
      B=U/MC2;
      NU=B+GA-1.;
      E=B/NU;
      KH=sqrt(B/(GA-1.));
      NS=1.;
      NB=exp(-4.*NS*KH*atan(1./KH));
      DB=1.-exp(-2.*PI*NS*KH);
      RECF=9614.*E*E*B/(GA-1.);
	  // Oleg AUNP=2.*KH/(1.-KH*KH-NU*NU);
//     SRK for dipole only (Bethe-Salpeter expression)
      SRK=RECF*NB/DB;
      AUSA=1.+3./4.*GA*(GA-2.)/(GA+1.)*(1.-log((1.+BETA)/(1.-BETA))/2./BETA/GA/GA);
//     Sauter cross. per electron
      SRKSA=3.77E-09*pow(Z,5)*BETA*GA/NU/NU/NU*AUSA/2.;

      if( SRK > SRKSA )  *SRKBEST=(SRK+SRKSA)/2.;
      else 		 *SRKBEST=SRKSA;

      *IRC = 0;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void reccorlm(double U, double *SRL, double *SRM, int *IRC ) {
//
//     REC subroutine for LM shell
//
      double   B, NU, E, LH, NS, NB, DB, RECF;
// **********************************************************************
//
//     start of LM REC calculation
      B=U/MC2;
      NU=B+GA-1.;
      E=B/NU;
      LH=sqrt(B/(GA-1.));
      NS=2.;
      NB=exp(-4.*NS*LH*atan(1./LH));
      DB=1.-exp(-2.*PI*NS*LH);
      RECF=9614.*E*E*B/(GA-1.);
//     dipole only
      *SRL=RECF*NB/DB*8.*(1.+6.*E+8.*E*E);
      *SRM=8./27.*(*SRL);
      *IRC = 0;
}

