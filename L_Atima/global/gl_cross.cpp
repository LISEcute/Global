#include <QtMath>
#include <stdio.h>
#include <stdlib.h>
#include "globallib.h"
#include "ion_potential.h"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
/*
C***********************************************************************
C     GLOBAL                                                           
C     A program to calculate and print out charge distributions        
C     of relativistic projectiles traversing solid or gaseous targets. 
C     A short description how to handle the program is given in the    
C     file global.readme. The physics basis of the program is described
C     in W.E. Meyerhof et al., Nucl. Inst. Meth. B (to be published)   
C                                                                      
C     W.E. Meyerhof, Stanford University, 1985 - 1996                  
C     E-mail: FE.WEM@FORSYTHE.STANFORD.EDU                             
C     Bertram Blank, CEN Bordeaux-Gradignan, 1993 - 1996               
C     Helmut Weick, GSI Helmholtzzentrum, h.weick@gsi.de  2015 - 2018  
C
C         ~2003 :   translated to cpp for LISE++, Oleg Tarasov, MSU
C          2015 :   ANSI-C version, H. Weick 
C     28.10.2016:   longitudinal ionization to zero for too low energies.
C                   otherwise leads to division by zero at some energy
C     31.10.2016:   init PP=0 in intpol, oterwise nonsense values
C     31.10.2016:   I_BRAN=0, may be good at some energies, but even wrong 
C                   direction at high energy (polarisation vs. screening), IONICRO
C     01.05.2017:   corrected interpolation, new lin. and log interpolation
C     20.06.2017:   made a subset without input and integration for cross sections
C                   at fixed energy, ignore double exchange and solve analytically
C                   -> equilib. charge-state distribution and mean values
C     06.07.2017:   extrapolated projectile ionization up to Z=120  
C***********************************************************************
*/

extern double   SKIB, BP1S,  P1S,  SKIS, BP1SS, P1SS,
SU2S, BP2SU, P2SU, SU2P, BP2PU, P2PU,
SS2S, BP2SS, P2SS, SS2P, BP2PS, P2PS,
SMIU, SMIS,  KIB,  KCB,  KIS, LIS,  MIS;

double BSA, GA, BETA, BSQ, D1a, D2, GAG, ETA;

void global_cross(int i_zf, int ne, int i_zm, double EN, double *d_s)
{
  double SC[4];

  double SRKBEST, SRLU, SRMU, LMCKB, KCBNR=0, KCBRA=0;
  double LCBNR, LCBRA, MCBNR, SRLS, SRMS, /*KCS,*/ LMCKS;
  double KCSNR, KCSRA, LIU=0, LCU=0, MCLU, LCS=0, /*MCLS,*/ MIU=0, MCU=0;
  double MCS=0, KRD, LRD, ZTFSQ, ZFALB2, ZMAX, ZAL, ZALB;
  double EXZAL, EXZALB, BMC, BMIU, BMIS;
  double KIBN, KISN, LIUN, LISN, MIUN, MISN, DKION, DLION /*, CSN*/;
  double LMCKBN, LMCKSN, LCUN, MCLUN, LCSN, /*MCLSN,*/ MCUN, MCSN ;
  double TCAPN, DLI, DMI, DLC, DMC;
  double ZF,ZT,ZP;
  double U1S, U1SS, U2S, U2P, U3S, U3D;

  double B[9] = { 0.,0., 0., 0., 0., .6, .3, 1., 1.};   // B[0:8]

  double KEX=2.,KMETA=100000., LEX=1.5, LMETA=10000.,LISSCR1=1.;
  double LISSCR2=4.,LIUSCR1=0.4,LIUSCR2=0.5,KISCR1=0.6,KISCR2=0.;
  double B3FAC=2.;
  int I, I_GAS, IRC, iCSoutput;
  int J = 28 ;
  double A[173];
  int I_BRAN = 1 ; // ALWAYS
  FILE *fp=NULL;
  const char *filename ;


  I_GAS = 0 ;  // gas or solid
  iCSoutput = 0 ; // =3 then to file
  IRC = 0 ;    // neg. means error

  // gas or solid
  /*      if((i_zm==1)||
         (i_zm==2)||
         (i_zm==7)||
         (i_zm==8)||
         (i_zm==9)||
         (i_zm==10)||
         (i_zm==17)||
         (i_zm==18)||
         (i_zm==36)||
         (i_zm==54)||
         (i_zm==86)
         ) I_GAS=1 ; */ // all as solids

  if(EN < 30) {
      printf("*** WARNING: E<30 MeV/u no x-sections, no charge exchange in matter\n");
      d_s[0]= 0. ;
      d_s[1]= 0. ;
      d_s[2]= 0. ;
      d_s[3]= 0. ;

      IRC = -1 ;
      return ;
    };

  ZF = (double) i_zf ;
  ZT = (double) i_zm ;

  U1S  = PRODAT[i_zf][1];
  U1SS = PRODAT[i_zf][2];
  U2S  = PRODAT[i_zf][3];
  U2P  = PRODAT[i_zf][4];
  U3S  = PRODAT[i_zf][5];
  U3D  = PRODAT[i_zf][6];
  

  // *** ioniz. cross. calc. at projectile energy
  //  ** data for ioniz. cross

  GA  = 1.+EN/931.494061;
  BSQ = 1.-1./GA/GA;
  BETA  = sqrt(BSQ);
  BSA=18780.*BSQ;

  ionicro(ZF, U1S, U1SS, U2S, U2P, U3S, U3D, &IRC);

  if( IRC != -8 && IRC != 0 && fp ) fprintf(fp,"\n<E>: RC from ionicro =%d\n",IRC);

  //
  // *** capt. cross calc. at projectile energies
  //  ** quantities used in non-radiative capture.
  D2=(GA-1.)/(GA+1.);
  D1a=sqrt(D2);
  GAG=GA/(GA+1.);
  ETA=ALPH/BETA;
  double ZTFAC=ZT*ZT+ZT;
  double DION=0.00231*ZT*ZT/ZF/BETA;
  int ISTATE ;

  for(ISTATE=1;  ISTATE<=6;  ISTATE++)    //L200

    //      calculate K,L,M ioniz. and capt. cross sections each unscreened
    //      and screened
    switch (ISTATE) {
      case 1 :                 // ***  capt. into K proj. states from K,L,M target states
        KIB = ZTFAC * SKIB;
        capcross(ZF, ZT, SC, &IRC);
        reccork(ZF,   U1S,  &SRKBEST, &IRC);
        reccorlm(U2S, &SRLU, &SRMU,   &IRC);
        KCB = SC[1] + ZT * SRKBEST;
        LMCKB = SC[2] + SC[3] + ZT * (SRLU + SRMU);
        if( I_BRAN == 1 ) KIB = KIB /pow(1.+ZT*BP1S,P1S);
        KCBNR = SC[1];
        KCBRA = ZT * SRKBEST;
        LCBNR = SC[2];
        LCBRA = ZT * SRLU;
        MCBNR = SC[3];
        break;

      case 2 :
        KIS = ZTFAC * SKIS;
        ZP = ZF-0.3;
        capcross(ZP, ZT, SC, &IRC);
        reccork(ZP, U1SS, &SRKBEST, &IRC);
        reccorlm(U2P, &SRLS, &SRMS, &IRC);
        // Oleg KCS = SC[1] + ZT * SRKBEST;
        LMCKS = SC[2] + SC[3] + ZT * (SRLS + SRMS);
        if( I_BRAN == 1 ) KIS = KIS /pow(1.+ZT*BP1SS,P1SS);
        KCSNR = SC[1];
        KCSRA = ZT * SRKBEST;
        break;

      case 3 :  // ***  capt. into L proj. states from K,L,M target states
        LIU = ZTFAC * (SU2S+3.*SU2P)/8.;
        ZP = ZF - 1.7;
        capcross(ZP, ZT, SC, &IRC);
        LCU = SC[2] + ZT * SRLU;
        MCLU = SC[3] + ZT * SRMU;
        if( I_BRAN == 1 )   LIU=ZTFAC*(SU2S/pow(1.+ZT*BP2SU,P2SU)+
                                       3.*SU2P/pow(1.+ZT*BP2PU,P2PU))/8.;
        break;

      case 4 :
        LIS = ZTFAC*(SS2S+3.*SS2P)/8.;
        ZP = ZF-3.;
        capcross(ZP, ZT, SC, &IRC);
        LCS = SC[2]+ZT*SRLS;
        //Oleg MCLS = SC[3]+ZT*SRMU;
        if( I_BRAN == 1 )     LIS=ZTFAC*(SS2S/pow(1.+ZT*BP2SS,P2SS)+
                                         3.*SS2P/pow(1.+ZT*BP2PS,P2PS))/8.;
        break;

      case 5 :   // ***  capt. into M proj. states from K,L,M target states
        MIU = ZTFAC*SMIU;
        ZP = ZF-8.;
        capcross(ZP, ZT, SC, &IRC);
        MCU = SC[3]+ZT*SRMU;
        break;

      case 6 :
        MIS = ZTFAC * SMIS;
        ZP = ZF - 19.;
        capcross(ZP, ZT, SC, &IRC);
        MCS = SC[3] + ZT * SRMS;
      };


  //L200
  //
  // *** calc. of B factors
  //
  if (ZF>2) {KRD=(LIU-KIB)/(KEX*LIU+KMETA);} else {KRD=0;}
  if (ZF>10){LRD=(MIU-LIU)/(LEX*MIU+LMETA);} else {LRD=0;}
  ZTFSQ = ZT*ZT/ZF/ZF;
  ZFALB2 =pow2(ZF*ALPH/BETA);
  B[0] = (1.+LRD)/(1+LISSCR1*ZTFSQ*(LISSCR2+ZFALB2));
  B[1] = (1.+LRD)/(1+LIUSCR1*ZTFSQ*(LIUSCR2+ZFALB2));
  B[2] = (1.+KRD)/(1+KISCR1 *ZTFSQ*(KISCR2 +ZFALB2));

  ZMAX = qMax(ZT,ZF);
  ZAL  = ZMAX*ALPH;
  ZALB = ZAL/BETA;
  EXZAL  = exp(ZAL);
  EXZALB = exp(ZALB);
  B[3] = (1+B3FAC*EXZALB) / (1+B3FAC*EXZAL);
  B[4] = 1.-KRD;
  BMC  = 1.-LRD;

  // do not exclude the few 4s states better use M-shell cross section instead
  //if( J < IMAX ) {} BMIU=BMIS=0;}
  //else {
  BMIU = B[1];
  BMIS = B[0] * 0.8;
  //     }

  if( I_GAS == 1 ) { //       gas target
      for(I=0; I<=8; I++) B[I] = 1.;
      BMIS = BMIU = BMC = 1.;
    }
  //
  // *** prepare final cross sections
  //
  KIBN = KIB*B[2];
  KISN = KIS*B[2];
  LIUN = LIU*B[1];
  LISN = LIS*B[0];
  MIUN = MIU*BMIU;
  MISN = MIS*BMIS;
  DKION = DION*B[5];
  DLION = DION*B[6];
  double KCBN = (KCBRA+KCBNR)*B[3];
  double KCSN = (KCSRA+KCSNR)*B[3];
  LMCKBN = LMCKB*B[4];
  LMCKSN = LMCKS*B[4];
  LCUN = LCU*B[4];
  MCLUN = MCLU*BMC;
  LCSN = LCS*B[4];
  // Oleg MCLSN = MCLS*BMC;
  MCUN = MCU*BMC;
  MCSN = MCS*BMC;
  double DCAPN = KCBNR/900000.*B[7];
  TCAPN = 1.3*DCAPN*DCAPN*B[8];
  // differences per electron between screened and unscreened
  DLI = (LISN-LIUN)/7.;
  DMI = (MISN-MIUN)/17.;
  DLC = (LCUN-LCSN)/7.;
  DMC = (MCUN-MCSN)/17.;

  // *** added for MOCADI charge-exchange straggling
  //     calculate cross secton for given charge state +/- 1 and +/- 2
  //     use values with solid state factor included

  // one electron ionization
  if(ne==0) // bare ion
    d_s[1] = 0. ;
  if(ne==1) // 1s, K-shell partially filled
    d_s[1] = KIBN ;
  if ((ne>=2)&&(ne<=10)) // 1s^2 + L-shell partially filled
    d_s[1] = 2*KISN + (ne-2)* ( (10-ne)/8.*LIUN + (ne-2)/8.*LISN );
  if ((ne>10)&&(ne<=28)) // K+L shell filled, M-shell partially
    d_s[1] = 2*KISN + 8*LISN + (ne-10)* ( (28-ne)/18.*MIUN + (ne-10)/18.*MISN ) ;
  if (ne>28) // even more bound electrons
    d_s[1] = 1.E9 ; // very large !

  // one electron capture
  if(ne==0) // capture bare
    d_s[2] = KCBN + LCUN + MCUN ;
  if(ne==1) // capture one electron in K-shell
    d_s[2] = 0.5*KCBN + LCUN + MCUN ;
  if(ne==2) // capture K-shell filled
    d_s[2] = LCUN + MCUN ;
  if ((ne>2)&&(ne<=10))  // capture K-shell filled, L-shell partially
    d_s[2] = (10-ne)/8.* ( (10-ne)/8.*LCUN +(ne-2)/8.*LCSN)+ MCUN ;
  if ((ne>10)&&(ne<=28)) // capture K,L-shell filled, M-shell partially
    d_s[2] = (28-ne)/18.*( (28-ne)/18.*MCUN +(ne-10)/18.*MCSN) ;
  if (ne>28)
    d_s[2] = 0. ;

  d_s[0] = 0. ;  // no double ionization
  d_s[3] = 0. ;  // no double capture

  //if(d_s[1]<0) {
  //  printf(" ioniz=%f, ne=%d KIB=%f KIS=%f %f %f \n",d_s[1],ne,KIB,KIS,SKIB,SKIS);
  //}
  //printf(" KCB=%f LCU=%f LCS=%f MCU=%f MCS=%f\n",KIBN,LCUN,LCSN,MCUN,MCSN);


  ////////////////////////////////////////////////////////////
  // prepare A() parameters - only needed for integration within GLOBAL


  A[0] = 0.;
  A[1] = 0.;
  A[2] = 0.;
  A[6] = 0.;
  A[7] = 0.;
  A[12] = 0.;
  A[3] = -(1.+DCAPN+TCAPN)*(KCBN+LMCKBN);
  A[4] = KIBN;
  A[10] =2.*KISN;

  for(I=2;  I<=9; I++)  A[6*I+4] = A[6*I-2]+LIUN+2*(I-2)*DLI;

  A[8]  = KCBN+LMCKBN;
  A[14] = KCSN/2.+LMCKSN;
  A[9]  = -A[4]-(1.+DCAPN+TCAPN)*A[14];

  for(I=3; I<=10; I++) A[6*I+2] = (LCUN-(I-3)*DLC)*(11-I)/8+MCLUN;

  if( J > 10 ) {
      for(I=10; I<=J-1; I++) A[6*I+4] = A[6*I-2]+MIUN+2.*(I-10)*DMI;
      for(I=11; I<=J  ; I++) A[6*I+2] = (MCUN-(I-11)*DMC)*(29-I)/18.;
    }

  for(I=0; I<=1; I++)   A[6*I+5] = DKION*A[6*I+10]*(I+1)/2.;
  for(I=2; I<=J-2; I++) A[6*I+5] = DLION*A[6*I+10]*(I+1)/2.;

  //     readjustment for L and M shells
  DCAPN /=2.;
  TCAPN /=2.;

  for(I=2; I<=J; I++)        A[6*I+1] = DCAPN*A[6*I-4];
  for(I=3; I<=J; I++)        A[6*I  ] = TCAPN*A[6*I-10];
  for(I=2; I<=J-3; I++)      A[6*I+3] = -A[6*I-2]-A[6*I-7]-A[6*I+8]*(1.+DCAPN+TCAPN);


  A[6*J-9] = -A[6*J-14]-A[6*J-19]-A[6*J-4]*(1.+DCAPN);
  A[6*J-3]=-A[6*J-8]-A[6*J-13]-A[6*J+2];
  A[6*J+3]=-A[6*J-2]-A[6*J-7];
  int IMIN = 6 * J+4;
  for(I=IMIN; I<=172; I++)  A[I]=0.;

  //

  // ***  list cross sections to file (debug option)

  if( iCSoutput == 3) {
      // if( 3 == 3) {
      filename = "global_cs-data.dat" ;
      fp=fopen(filename,"wt");
      if(!fp) {
          IRC=-1;
          return ;
        }

      fprintf(fp,"\n<I>: EN = %.2f MeV/u\n",EN);

      fprintf(fp,"\n KIB       LIU       MIU       KCB       LCU       MCU       DbleKion  LMCKB     DbleLion");
      fprintf(fp,"\n KIS       LIS       MIS       KCBNR     LCBNR     MCBNR     DbleCapt  TOTCAP    TrplCapt");
      fprintf(fp,"\n                               KCBRA     LCBRA     MCLU ");

      fprintf(fp,"\n %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e",
              KIBN,LIUN,MIUN,KCBN,LCUN,MCUN,DKION,LMCKBN,DLION);

      fprintf(fp,"\n %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e %-9.2e ",
              KISN,LISN,MISN,KCBNR*B[3],LCBNR*B[4],MCBNR*BMC,DCAPN*2.,KCBN+LMCKBN,TCAPN*2.);

      fprintf(fp,"\n                               %-9.2e %-9.2e %-9.2e",
              KCBRA*B[3],LCBRA*B[4],MCLUN);


      fprintf(fp,"\n Factors:");
      fprintf(fp,"\n %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f",
              B[2], B[1], BMIU, B[3], B[4], BMC, B[5], B[4],  B[6]);
      fprintf(fp,"\n %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f    %-6.4f              %-6.4f",
              B[2], B[0], BMIS, B[3], B[4], BMC, B[7], B[8]);

      fprintf(fp,"\n                               %-6.4f    %-6.4f    %-6.4f  \n",
              B[3], B[4],BMC);

      fclose(fp);
    }

  /*

 KIB    ioniz. cross sec. for one 1s        electron
 KIS      "     "     "    "  one (1s)^2       "
 LIU      "     "     "    "  one 2s or p      "
 LIS      "     "     "    "  one (2s or p)^6  "
  the L-shell cross section is an average of s and p states.
 MIU      "     "     "    "  one 3s,p,d      electron
 MIS      "     "     "    "  one (3s,p,d)^18    "
  the M-shell cross section is an average of s, p, and d states.

 KCB    NR+RAD capt. cross into empty proj. K shell         from target KLM
 LCU                   "                    L      , but full K shell    "
 MCU                   "                    M      ,  "       KL shells  "
 KCBNR  NR part of KCB
 LCBNR      "      LCU
 MCBNR      "      MCU
 KCBRA  RAD part of KCB
 LCBRA      "       LCU
 LMCKB  NR+RAD capt into empty LM shells if K shell is empty
 MCLU              "            M           L
 TOTCAP=KCB+LMCKB for the total capture cross section.

*/

}

