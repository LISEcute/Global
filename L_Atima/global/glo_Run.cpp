#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctime>
#include <iostream>

//#include "local_constant.h"
#include "globallib.h"
#include "glo_targetData.h"
#include "ion_potential.h"
//#include "gl_data_cst.h"

#include <QtMath>

//#include <QDebug>
#define SizeBuf 200


extern void ELEMENT(double &Z, double &A, char *CZ, int IOPT, int &IRC );
extern double RANGE(double A, double Z, double AT, double ZT, double T);
extern double ENERGY(double A, double Z, double AT, double ZT, double RG);
extern char *GetDoubleFromString(double &V, char *s);


extern void CROSS(double ZF, double ZT, int iCSoutput,
			double U1S, double U1SS, double U2S, double U2P, double U3S,
			double U3D, double EN,   int J, int I_GAS, int &IRC, FILE *fp);



extern double ar_fk0[];
extern double ar_fl1[];
extern double ar_fl2[];
extern double ar_fm1[];
extern double ar_fm2[];
extern double ar_fm3[];

double A[173];     // it is ok  A[0:172]

extern double   SKIB, BP1S,  P1S,  SKIS, BP1SS, P1SS,
SU2S, BP2SU, P2SU, SU2P, BP2PU, P2PU,
SS2S, BP2SS, P2SS, SS2P, BP2PS, P2PS,
SMIU, SMIS,  KIB,  KCB,  LIS,  MIS, KIS;


int RunGlobalLocal(double *Obmen, char *filename, bool option_read_data_file,
			 int fast, char *Global_version,
			 double gAF, double gZF, double gAT, double gZT,
			 double gDTARGET, double gEN0, int  gQIN,
			 int   I_WR, int   iOption, int   iCSoutput,
			 int   iLoop, int N_Steps, int  Qshow,
			 int    DELZF, double DELE,
			 int    DELQ, int DELZT, double DELDT);

int FillArrayFromFile(const char *fname, double *ar);

int GlobalExtern(double *qtab, double *Inform, int fast, double gAF, double gZF, double gAT,
		     double gZT, double gDTARGET, double gEN0, int gQIN, int iOption);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//     MF version 01/97  (original version.  no modifications)
//*****************************************************************
//     f77 global.f -L/usr/local/cern/new/lib -lkernlib -lmathlib
//                  -Nl100 -o global
//*****************************************************************
// HOUSE KEEPING REMARKS
//
//     Version 08/96: includes EN adjustment, changes in Output,
//                    new I_WRI fro WRITE frequency
//             09/96: contains D_EQ criterium to check equilibrium
//                    thickness
//                    DELT calculated from LIS or MIS, if too big
//                    devided by 10
//                    DELT no longer in MENUE
//                    Energy criterium for new x-section calculation
//                        log10(en)**2
//                    no solid-state factors for gas targets
//             10/96: J=max+5 (before +2)
//                    all variables defined
//                    new subroutine uppercase
//                    "cleaned" version
//                    modified SUBs PILO, PILO10, PRSO
//             11/96: modified output: Eout, D_eq in colums
//                    int-step / 10
//             01/97: debugged from test persons, mod. out-freq.,
//                    correct A,Z,Symbol, minor buggs corrected
//             02/97: I_WRITE_FAC added, but should probably be removed ;
//                    DELT removed from MENUE subroutine
//                    charge changed to Z for the loop question
//                    "Target step size too large" is written also on file;
//**********************************************************************
//     GLOBAL
//     A program to calculate and print out charge distributions
//     of relativistic projectiles traversing solid or gaseous targets.
//     A short description how to handle the program is given in the
//     file global.readme. The physics basis of the program is bed*;
//     in W.E. Meyerhof et al., Nucl. Inst. Meth. B (to be published)
//
//     W.E. Meyerhof, Stanford University, 1985 - 1996
//     E-mail: FE.WEM@FORSYTHE.STANFORD.EDU
//     Bertram Blank, CEN Bordeaux-Gradignan, 1993 - 1996
//     E-mail: blank@cenbg.in2p3.fr
//**********************************************************************
//     Notation: ZF,AF:  proj. Z,A;  ZT,AT:  target Z,A
//               X: calculated charge fractions
//               DX: change of charge-state fraction in integration ;
//               Exx: incident-energy-related variables (almost !)*;
//               Dxx: target-thickness-related variables
//               Qxx: charge-states-related variables
//               TXx: integration step sizes
//               Uxx: ionisation potentials for different shells
//**********************************************************************
//

int RunGlobalLocal(double *Obmen, char *filename, bool option_read_data_file,
                   int fast, char *Global_version,
                   double gAF, double gZF, double gAT, double gZT,
                   double gDTARGET, double gEN0, int  gQIN,
                   int   I_WR, int   iOption, int   iCSoutput,
                   int   iLoop, int N_Steps, int  Qshow,
                   int    DELZF, double DELE,
                   int    DELQ, int DELZT, double DELDT)

{

  //=========================================  LOCAL
  double AF,ZF,AT,ZT,NT,DE,DTARGET,EN0, EN;
  double AF0,ZF0,AT0,ZT0,DTARGET0;
  int   QIN, QIN0, I_EN, I_WRI, kkl;

  double ZFMIN, EMIN, QMIN, ZTMIN, DTMIN, ENMAX;
  int JQST[IMAX];
  //======================================
  const char *message_low_energy="\n<E>: Energy lower than 30 MeV/u!\nCalculations stopped!";


  FILE *fp=nullptr;

  if(filename) {
      fp=fopen(filename,"wt");
      if(!fp) return -697;
    }

  const char* fk0 = "fk0";
  const char* fl1 = "fl1";
  const char* fl2 = "fl2";
  const char* fm1 = "fm1";
  const char* fm2 = "fm2";
  const char* fm3 = "fm3";

  if(option_read_data_file) {
      if(FillArrayFromFile(fk0,ar_fk0)<0) return -2;
      if(FillArrayFromFile(fl1,ar_fl1)<0) return -3;
      if(FillArrayFromFile(fl2,ar_fl2)<0) return -4;
      if(FillArrayFromFile(fm1,ar_fm1)<0) return -5;
      if(FillArrayFromFile(fm2,ar_fm2)<0) return -6;
      if(FillArrayFromFile(fm3,ar_fm3)<0) return -7;
    }

  fast=qMin(2,qMax(fast,0));

  int Result=1;
  if(iCSoutput!=3) iCSoutput=2;
  int loopLocal=0;       // for LOCAL use
  int I_GAS;
  int IRC, IP, IT, J, I, IAA, NumberCickles;
  bool equilibriumReached=false;
  char CQST[IMAX][15];
  char strline[81]="-------------------------------------------------------------------------------\n";
  double DTARG=0, QMEAN=0, EOUT;
  double ENKEEP, AUX, DELT=0, D_EQ=0, RG,EDIFF, SX,/*EA,ZSQ,*/ dQ=0,t;
  double E_final, E_init, DTARGET_print;
  //double D_Eq_final;
  double TX, TXA[6];  // it is ok
  long unsigned  int LLOW[6], LHIGH[6], LINV[6];
  double X[31], DX[31], SA[29], FA[29];   //it is ok X[0:30], DX[0:30], SA[0:28], FA[0:28]
  double coef_TXA[6]={1,10,100,1000,10000,100000};
  char CPRO[10];
  char CTAR[10];

  long unsigned  int  total_loop=0, total_cross=0;

  int I_WRITE_FAC=0;
  bool flagToStop, bool_flag2, loop_for_opt_EqEnd;
  long unsigned int imaxloop, IR=0;
  int Qstop=qMin(int(gZF),qMin(Qshow+9,IMAX-1));
  int Len;
  double U1S,U1SS,U2S,U2P,U3S,U3D;

  // struct date mdate;   struct time mtime;
  // getdate(&mdate);     gettime(&mtime);
  time_t ti = time(0);
  struct tm * now = localtime(&ti);



  /*
C     IF( I_BRAN .EQ. 0 ) CBRANDT = 'N'
C     IF( I_BRAN .EQ. 1 ) CBRANDT = 'Y'
C     CALL PYES(
C    # 'Use binding-polarization factor for x-sections (y/n)', CBRANDT)
C     IF( CBRANDT .NE. 'N' .AND. CBRANDT .NE. 'n' ) I_BRAN = 1
*/
  double QMEANMX=0;
  bool optionFlagEquil  = (iOption==opt_EqIni   || iOption==opt_EqEnd);
  bool optionFlagNorm   = (iOption==opt_NormIni || iOption==opt_NormEnd);
  bool optionFlagEfinal = (iOption==opt_EqEnd   || iOption==opt_NormEnd);
  bool FlagTargetPrint=false;

  AF  = gAF;
  AT  = gAT;
  ZF  = ZFMIN = gZF;
  EN0 = EMIN  = gEN0;

  if(gQIN==777) {QIN = QMIN  = 0.;}
  else          {QIN = QMIN  = gQIN;}

  ZT  = ZTMIN = gZT;
  DTARGET = DTMIN = gDTARGET;

  if( optionFlagEquil && iLoop==5) iLoop=0;

  for(int i=0; i<IMAX; i++) JQST[i]=i;
  //=================================================== start

  DTARGET0 = DTARGET+N_Steps*DELDT;
  ENMAX = EN0+N_Steps*DELE;
  ZF0 = ZF;
  AF0 = AF;
  QIN0 = QIN;
  ZT0 = ZT;
  AT0 = AT;
  loopLocal = 0;
  if( optionFlagEquil ) gDTARGET = 100000.;

  I_WRI = pow_int(10,I_WR+1);
  //=================================================================
  //=================================================================
L5:

  EN = EN0 = gEN0;       /// init for each loopLocal  -> EN , EN0  can be changed
  DTARGET = gDTARGET;    //  init for each loopLocal  -> DTARGET  can be changed

  I_EN = 0;


  switch(iLoop)
    {
    case  1  :
      ZF = ZFMIN + double(loopLocal * DELZF);
      ELEMENT( ZF, AF, CPRO, 3, IRC );
      break;
    case  2 :
      EN  = EMIN + double(loopLocal * DELE);
      EN0 = EN;
      break;
    case  3 :
      QIN = int( QMIN + double(loopLocal * DELQ));
      //              J = 28;
      break;
    case  4 :
      ZT = ZTMIN + double(loopLocal * DELZT);
      ELEMENT( ZT, AT, CTAR, 1, IRC );
      break;

    case  5 :
      DTARGET = DTMIN + double(loopLocal * DELDT);
      break;
    };

  IP = ZF;
  //J = qMin(28,max(IMAX,qMin(9,QIN+5)));   //OLEG
  J = qMax(1,qMin(int(ZF)-1, qMax(IMAX,qMin(9,QIN+5))));

  IT = int(ZT);
  if( ZT == 6.6 ) IT = 97;

  if( loopLocal == 0 && fp )
    {
      fputs(strline,fp);
      fputs("       GLOBAL: Q-states of heavy ions behind matter layers\n",fp);
      fputs(strline,fp);
    }

  if( ( iCSoutput == 1 || iCSoutput == 2 ) && loopLocal == 0  && fp )
    {
      fprintf(fp,"***** Global *****    %s    ****       %02d-%02d-%04d   %02d:%02d:%02d\n",
              Global_version, now->tm_mday, now->tm_mon+1, now->tm_year+1900,
              now->tm_hour, now->tm_min, now->tm_sec);
      fputs(strline,fp);
    }

  //
  //=============================================================
  // *** Final charge fraction calc. and printout
  //
  for( I=Qshow; I<=Qstop; I++)
    if( JQST[I] >= 0 )      sprintf(&CQST[I][0],"   Q(%-2d)  ",JQST[I]);
    else                    sprintf(&CQST[I][0]," ");

  //==============================================
  // =============================================   loopLocal==0
  //
  const char *L31 =" (Z=%d A=%.0f Qe=%d) at E=%.1f MeV/u   on (Z=%d  A=%.1f  D=%.3g mg/cm^2)";
  const char *L322=" (Projectile: Qe=%d) at E=%.1f MeV/u   on (Z=%d  A=%.1f";
  const char *L333=" (Z=%d A=%.0f Qe=%d)  on (Z=%d  A=%.1f";
  const char *L344=" (Z=%d A=%.0f) at E=%.1f MeV/u  on (Z=%d  A=%.1f";
  const char *L355=" (Z=%d A=%.0f Qe=%d) at E=%.1f MeV/u  on (";
  const char *L36 =" (Z=%d A=%.0f Qe=%d) at E=%.1f MeV/u  on (Z=%d  A=%.1f";
  const char *L3t ="  D=%.4g mg/cm^2)";


  if( loopLocal == 0 && fp )
    {
      DTARG = ( optionFlagEquil ? 0 : DTARGET);
      if( iCSoutput == 2 )
        {
          fprintf(fp,"\n");

          switch(iLoop)
            {
            case 0 :
            case 5 : fprintf(fp,L36, int(ZF), AF, QIN, EN0, int(ZT), AT); break;

            case 1 : fprintf(fp,L322,        QIN, EN0, int(ZT), AT); break;
            case 2 : fprintf(fp,L333,int(ZF), AF, QIN, int(ZT), AT);  break;
            case 3 : fprintf(fp,L344,int(ZF), AF, EN0, int(ZT), AT); break;
            case 4 : fprintf(fp,L355,int(ZF), AF, QIN, EN0); break;
            };

          if( iLoop != 5 && !optionFlagEquil) fprintf(fp,L3t,DTARG);
          else                            fprintf(fp,")");
        }
    }


  if( iCSoutput == 3 && fp )
    {
      fprintf(fp,"\n");
      fprintf(fp,L31, int(ZF),AF, QIN, EN0,int(ZT),AT, DTARG);
      fprintf(fp,"\n");
    };


  // ============================================   loopLocal==0

  if( loopLocal == 0  && fp)
    {
      fprintf(fp,"\n\n");
      if(  (optionFlagNorm || optionFlagEquil ||  ( iOption == opt_Evolution && iLoop == 0 ))  &&
           iCSoutput!=3 ) {


          const char *charge_str[] = {"Q-states at target exit (E init):",
                                      "Q-states at target exit (E final):",
                                      "Equilibrium Q-states (E init):",
                                      "Equilibrium Q-states (E final):",
                                      "Q-states evolution:"};



          if( iCSoutput == 0 || iCSoutput == 2 ) fputs(charge_str[iOption],fp);
          fprintf(fp,"\n\n");


          int flag_print_charge=0;


          if( iLoop == 0 && iOption == opt_Evolution ) flag_print_charge=1;
          else if( iLoop == 0 && optionFlagEquil )               flag_print_charge=2;
          else if( iLoop == 0 )                flag_print_charge=3;
          else if( iLoop == 1 )                flag_print_charge=4;
          else if( iLoop == 2 )                flag_print_charge=5;
          else if( iLoop == 3 )                flag_print_charge=6;
          else if( iLoop == 4 )                flag_print_charge=7;
          else if( iLoop == 5 )                flag_print_charge=8;


          const char *char_energy[2]= {"out    ","init   "};

          const char *strff[] = {
            " ",                                    // 0
            "D(mg/cm2)  ",                          // 1
            "D_eq(mg/cm2) E",                // 2
            "D(mg/cm2) D_eq    E",           // 3
            "  Z_p     D_eq    E",           // 4
            " E(MeV/u) D_eq    E",           // 5
            "  Q       D_eq    E",           // 6
            "  Z_t     D_eq    E",           // 7
            "D(mg/cm2) E"                    // 8
            //                "                              ",       // 9
            //                "D_eq(mg/cm2) Einit   "                 // 10
          };

          Len = (int)strlen(strff[flag_print_charge]);
          if(iOption != opt_Evolution) Len+=7;     // for char_energy

          for(I=0; I<Len; I++) fprintf(fp,"-");
          fprintf(fp,"----------------");
          for(I=Qshow; I<=Qstop; I++) fprintf(fp,"-------");
          //----------------------------------------------
          fprintf(fp,"\n%s",strff[flag_print_charge]);

          if(iOption != opt_Evolution)
            fputs(char_energy[optionFlagEfinal? 1 : 0],fp);


          fprintf(fp,"Qmean  dQ    *| ");
          for(I=Qshow; I<=Qstop; I++) fprintf(fp,"%s ",CQST[I]);
          //----------------------------------------------
          fprintf(fp,"\n");
          for(I=0; I<Len; I++) fprintf(fp,"-");
          fprintf(fp,"-------------*|-");
          for(I=Qshow; I<=Qstop; I++) fprintf(fp,"-------");
        }
    }


  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ///  for Energy Final and NORMAL mode
  DTARGET_print = DTARGET;
  E_init = gEN0;

  if(iOption == opt_NormEnd)
    {
      E_final = gEN0;
      RG   = RANGE( AF, ZF, AT, ZT, gEN0);
      RG  += DTARGET;
      EN0 = EN = ENERGY( AF, ZF, AT, ZT, RG);
      E_init = EN0;
    }

  if(iOption == opt_NormIni)
    {
      RG   = RANGE( AF, ZF, AT, ZT, gEN0);
      RG  -= DTARGET;
      E_final = ENERGY( AF, ZF, AT, ZT, RG);
      if( E_final <= 30. )  {
          I_EN = 1;
          if(fp) fputs (message_low_energy, fp);
          goto L9000;
        }
    }


  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


  //========================================================

  if( iOption == opt_Evolution && iLoop != 0 )
    {
      if( iCSoutput == 2  &&  fp)
        {
          if( loopLocal == 0 ) fprintf(fp,"Q-state evolution:");

          fprintf(fp,"\n\n");

          switch(iLoop) {
            case 1 : fprintf(fp," Projectile: Z = %d, A = %d",int(ZF+0.5), int(AF+0.5)); break;
            case 2 : fprintf(fp," Energy:  E = %.1f MeV/u",EN); break;
            case 3 : fprintf(fp," Q-state: Q = %d+",QIN); break;
            case 4 : fprintf(fp," Target: Z = %d, A = %d",int(ZT+0.5), int(AT+0.5)); break;
            case 5 : fprintf(fp," Thickness: D = %.4g mg/cm^2",DTARGET_print); break;
            };

          fprintf(fp,"\n------------------------------");
          for(I=Qshow; I<=Qstop; I++) fprintf(fp,"-------");


          fprintf(fp,"\nD(mg/cm^2)    Qmean  dQ    *| ");
          for(I=Qshow; I<=Qstop; I++) fprintf(fp,"%s ",CQST[I]);

          fprintf(fp,"\n---------------------------*|-");
          for(I=Qshow; I<=Qstop; I++) fprintf(fp,"-------");
        }
    }

  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  //WWWWWWWWWWWWWWWWWWWWW      start           WWWWWWWWWWWWWWWWWWWWW
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  ZF   = PRODAT[IP][0];

  U1S  = PRODAT[IP][1];
  U1SS = PRODAT[IP][2];
  U2S  = PRODAT[IP][3];
  U2P  = PRODAT[IP][4];
  U3S  = PRODAT[IP][5];
  U3D  = PRODAT[IP][6];

  ZT   = TARDAT[IT][0];
  NT   = TARDAT[IT][1];
  DE   = TARDAT[IT][2];


  if( NT < 0.01 ) I_GAS = 1;
  else            I_GAS = 0;
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  // *** Determine equilibrium thickness and step size
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  //
  loop_for_opt_EqEnd=true;

label_loop_for_opt_EqEnd:

  ENKEEP = EN;

  CROSS(ZF, ZT, iCSoutput, U1S, U1SS, U2S, U2P, U3S, U3D, EN, J, I_GAS, IRC, fp);


  if( IRC != 0 ) goto L9000;
  D_EQ = DE * 100. / NT * 4.6 / (KIB+ KCB/2.);

  ///////////////////////////////////   OLEG mode   opt_EqEnd
  if(iOption == opt_EqEnd)
    {
      RG   = RANGE( AF, ZF, AT, ZT, gEN0);
      RG  += D_EQ;
      EN   = ENERGY( AF, ZF, AT, ZT, RG);
      if(loop_for_opt_EqEnd) {
          loop_for_opt_EqEnd=false;
          goto label_loop_for_opt_EqEnd;
        }
      EN0=EN;
      E_init = EN0;
    }
  ////////////////////////////////////////////////

  //////////////////////////////////////   Oleg mode
  //////////////////////////////////////   if thickness is more larger than D_EQ

  /*              comment on 10/30/2012

if( 3.*D_EQ < DTARGET && optionFlagNorm) {

      // for final energy
      ENKEEP = EN;
      CROSS(ZF, ZT, iCSoutput, U1S, U1SS, U2S, U2P, U3S, U3D, E_final, J, I_GAS, IRC, fp);
      if( IRC != 0 ) goto L9000;
      D_EQ = DE * 100. / NT * 4.6 / (KIB+ KCB/2.);
      //-----------------------------------
      RG   = RANGE( AF, ZF, AT, ZT, E_final);
      RG  += 3.*D_EQ;
      DTARGET = 3.*D_EQ;
      FlagTargetPrint = true;
      EN0= EN = ENERGY( AF, ZF, AT, ZT, RG);
      }
   */
  ///////////////////////////////////////////

        if( J <  10 ) AUX = DE / NT / KIS;
  else  if( J <  28 ) AUX = DE / NT / LIS;
  else                AUX = DE / NT / MIS;



  /*   OLEG
      if( J <= IMAX ) AUX = DE / NT / LIS;
      else            AUX = DE / NT / MIS;

*/

  if (AUX < 0.0003) AUX = 0.0001;   // Oleg
  else if (AUX < 0.003 ) AUX = 0.001;   // Oleg
  else if (AUX < 0.03  ) AUX = 0.01;
  //      else if (AUX < 0.08  ) AUX = 0.05;
  else if (AUX < 0.30  ) AUX = 0.1;
  //      else if (AUX < 0.80  ) AUX = 0.5;
  else                   AUX = 1.;

  DELT = 0.01 * AUX;
  if ( NT < 0.01 ) DELT = 0.002 * AUX;

  DELT/=10;    //OLEG

  if(fast>0 || gDTARGET>200) DELT*=10;   // v.4.6.20
  if(fast>1 || gDTARGET>2000) DELT*=10;  // v.4.6.20


  // *** RESET
L140:

  for( I=0; I<=30; I++) {  X[I] =DX[I]=0;}
  for( I=0; I<=28; I++) {  SA[I]=FA[I]=0;}

  if(gQIN==777)
    for(I=0; I<28; I++) X[I] = Obmen[I];  // OT  -- inital charge states
  else
    X[QIN]=1.;


  I_WRITE_FAC = 1;
  //
  // *** start of charge-fraction integration
  //     RESTART FOR NEGATIVE CHARGE STATES

  TX=NT*DELT/100./DE;

  for(int I=0; I<6; I++)
    {
      TXA[I] = TX*coef_TXA[I];
      LINV[I]= coef_TXA[I];
    }

  kkl=1;
  NumberCickles=100;
  //

  imaxloop = int(qMin((unsigned long)0xFFFFFFFF,(unsigned long)(DTARGET/DELT*2.1)));

  LLOW[0] = 1;                        LHIGH[0]=kkl*NumberCickles;  NumberCickles*=10;
  LLOW[1] = LHIGH[0]+coef_TXA[1];     LHIGH[1]=kkl*NumberCickles;  NumberCickles*=10;
  LLOW[2] = LHIGH[1]+coef_TXA[2];     LHIGH[2]=kkl*NumberCickles;  NumberCickles*=10;
  LLOW[3] = LHIGH[2]+coef_TXA[3];     LHIGH[3]=kkl*NumberCickles;  NumberCickles*=10;
  LLOW[4] = LHIGH[3]+coef_TXA[4];     LHIGH[4]=kkl*NumberCickles;  NumberCickles*=10;
  LLOW[5] = LHIGH[4]+coef_TXA[5];     LHIGH[5]=imaxloop;

  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  GLOBAL LOOP begin

  for(IAA=0; IAA<6; IAA++)
    {
      TX=TXA[IAA];

      for( I=0;  I<=J;  I++)
        {
          SA[I]=A[6*I]+A[6*I+1]+A[6*I+2]+A[6*I+3]+A[6*I+4]+A[6*I+5];
          SA[I]*=TX;
          FA[I]=1.+SA[I]/2.+SA[I]*SA[I]/6.+SA[I]*SA[I]*SA[I]/24.;
        }


      //===============================================================  L205 loop begin
      //
      for(IR=LLOW[IAA];  IR<=LHIGH[IAA];  IR+=LINV[IAA])  //L205
          {
          //    Adjust projectile energy, if necessary
          total_loop++;
          if(IR >= imaxloop)
            if(fp) fprintf(fp,"\n<I> Max.step number %lu reached!\n", IR);

          RG  = RANGE( AF, ZF, AT, ZT, EN0);
          RG -= DELT*IR;
          EN  = ENERGY( AF, ZF, AT, ZT, RG);

          if( EN <= 30. )
            {
              I_EN = 1;
              if(fp) fputs(message_low_energy,fp);
              goto L9000;
            }

          EDIFF = 1.- pow2(log10(EN))/ 200.;
          if( EN < ENKEEP * EDIFF || ENKEEP < EN )
            {
              ENKEEP = EN;
              total_cross++;
              CROSS(ZF, ZT, iCSoutput, U1S, U1SS, U2S, U2P, U3S, U3D, EN, J, I_GAS, IRC, fp);
              if( IRC != 0 )
                goto L9000;

              for( I=0;  I<=J;  I++)
                {
                  SA[I] = A[6*I]+A[6*I+1]+A[6*I+2]+A[6*I+3]+A[6*I+4]+A[6*I+5];
                  SA[I]*= TX;
                  FA[I] = 1.+SA[I]/2.+SA[I]*SA[I]/6.+SA[I]*SA[I]*SA[I]/24.;
                }

              D_EQ = DE * 100. / NT * 4.6 / (KIB+ KCB/2.);
            }

          //        EA=0.;
          SX=0.;
          //        ZSQ=0.;
          DX[0] = TX*(A[3 ]*X[0] + A[4 ]*X[1] + A[5 ]*X[2]                        )*FA[0];
          DX[1] = TX*(A[8 ]*X[0] + A[9 ]*X[1] + A[10]*X[2] + A[11]*X[3]           )*FA[1];
          DX[2] = TX*(A[13]*X[0] + A[14]*X[1] + A[15]*X[2] + A[16]*X[3]+A[17]*X[4])*FA[2];

          for( I=3;  I<=J;  I++)
            DX[I]=TX*(A[6*I  ]*X[I-3]+
                A[6*I+1]*X[I-2]+
                A[6*I+2]*X[I-1]+
                A[6*I+3]*X[I  ]+
                A[6*I+4]*X[I+1]+
                A[6*I+5]*X[I+2]
                )*FA[I];


          for( I=0;  I<=J;  I++)  if( fabs(DX[I]) < 9.999999E-21 ) DX[I] = 0.;

          for( I=0;  I<=J;  I++)
            {
              X[I]+=DX[I]; if( fabs(X[I]) < 1.E-20 ) X[I] = 0.;
              SX  +=X[I];  if( fabs(SX)   < 1.E-20 ) SX = 0.;
            }

          for( I=0;  I<=J;  I++)
            {
              if(SX!=0) X[I]/=SX;
              if( fabs(X[I]) < 1.E-20 ) X[I] = 0.;
              /*	     EA+=I*X[I];
             ZSQ+=(ZF-I)*(ZF-I)*X[I];*/
            }

          //   CHECK CHARGE STATES
          for( I=0; I<=J; I++)
            if( X[I] < -1.E-5)
              {
                DELT /= 10.;
                I_WRITE_FAC *= 10;
                ENKEEP = 0.;
                goto L140;
              }

          //  **   CHECK EQUILIBRIUM CHARGE STATES
          equilibriumReached = (DELT*IR >= D_EQ);


          QMEAN = 0;
          for(I=0; I<=J; I++) QMEAN += X[I]*(ZF-I);

          dQ=0;
          for(I=0; I<=J; I++)
            {
              t = ZF-I-QMEAN;
              dQ += X[I]*t*t;
            }
          dQ = (dQ>0 ? sqrt(dQ)  : 0);

          //
          //  **   OUTPUT

          if( iCSoutput == 2 && fp)
            {
              if( (optionFlagNorm && DELT*IR >= DTARGET ) ||
                  (optionFlagEquil && equilibriumReached )
                  )
                {

                  if(optionFlagNorm && iLoop==5)
                    {
                      if ( QMEAN >= QMEANMX ) QMEANMX = QMEAN;
                    }

                  flagToStop  = (iLoop==0 && optionFlagEquil && equilibriumReached);
                  bool_flag2 = (iLoop==0 && optionFlagNorm && FlagTargetPrint);

                  if( flagToStop   )  fprintf(fp,"\n%-9.2f  %-8.2f ",DELT*IR, (optionFlagEfinal ? E_init : EN));
                  else if( bool_flag2  )  fprintf(fp,"\n%-9.2f",DTARGET_print);
                  else if( iLoop == 0 )  fprintf(fp,"\n%-9.2f",DELT*IR);
                  else if( iLoop == 1 )  fprintf(fp,"\n   %-3.0f   ",ZF);
                  else if( iLoop == 2 )  fprintf(fp,"\n  %-6.1f ",EN0);
                  else if( iLoop == 3 )  fprintf(fp,"\n  %-3d    ",QIN);
                  else if( iLoop == 4 )  fprintf(fp,"\n   %-3.0f   ",ZT);
                  else if( iLoop == 5 )  fprintf(fp,"\n%-9.2f",DTARGET);

                  if(!flagToStop && iLoop != 5) fprintf(fp," %-7.2f %-6.1f ",D_EQ, (optionFlagEfinal ? E_init :EN ));

                  if( iLoop == 5 ) fprintf(fp," %-6.1f ", (optionFlagEfinal ? E_init : EN));

                  fprintf(fp," %-6.2f %-5.2f *|  ",QMEAN, dQ);

                  for(I=Qshow; I<=Qstop; I++)
                    if(X[JQST[I]] >1e-15)   fprintf(fp,"%9.3e  ",(X[JQST[I]]));
                    else                   fprintf(fp,"    0      ");

                  //                                if(X[JQST[I]] >1e-5)   fprintf(fp,".%04d  ",int(X[JQST[I]]*10000));
                  //                                else                   fprintf(fp,".0     ");

                }

              else  if( iOption == opt_Evolution  )
                {
                  if( ((IR)%( I_WRI*I_WRITE_FAC )) == 0 || DELT*IR >= DTARGET ) {

                      fprintf(fp,"\n%-9.2f ",DELT*IR);


                      if( iLoop != 0 ) fprintf(fp,"   ");
                      fprintf(fp," %-6.2f %-5.2f *|  ",QMEAN, dQ);

                      for(I=Qshow; I<=Qstop; I++)
                        if(X[JQST[I]] >1e-15)   fprintf(fp,"%9.3e  ",(X[JQST[I]]));
                        else                   fprintf(fp,"    0      ");

                      //                                if(X[JQST[I]] >1e-5)   fprintf(fp,".%04d  ",int(X[JQST[I]]*10000));
                      //                                else                   fprintf(fp,".0     ");

                    }
                }
            }


          if( (optionFlagNorm && DELT*IR >= DTARGET ) ||
              (optionFlagEquil && equilibriumReached) ||
              (iOption == opt_Evolution && DELT*IR >= DTARGET ))
            goto L9000;

        }
      //===============================================================  L205 loop end
    }
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  GLOBAL LOOP end
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
L9000:


  //
  if( !equilibriumReached && iLoop > 0 && optionFlagEquil )
    if(iCSoutput == 2 && fp) fprintf(fp,"\n            <W>: Equilibrium target thickness not reached!\n");

  //
  RG = RANGE( AF, ZF, AT, ZT, E_init);
  RG  -= IR * DELT;
  EOUT = ENERGY( AF, ZF, AT, ZT, RG);
  if(EOUT<30)I_EN=1;

  flagToStop=false;          // stop

  if( iLoop != 0 && I_EN==0)
    {
      flagToStop=true;          // continue
      loopLocal++;

      switch(iLoop)
        {

        case 1 :   if( ZF+DELZF > 96.) {
              ZF = ZF0;
              AF = AF0;
              flagToStop=false; // stop
            }
          break;

        case 2 :   if(EMIN + double((loopLocal+1) * DELE) > ENMAX) {
              EN = ENMAX;
              EN0 = ENMAX;
              flagToStop=false; // stop
            }
          break;

        case 3 :   if( QIN+DELQ >= int(ZF/2.) || QIN+DELQ > 28. ) {
              QIN = QIN0;
              flagToStop=false; // stop
            }
          break;

        case 4 :   if( ZT+DELZT > 96. ) {
              ZT = ZT0;
              AT = AT0;
              flagToStop=false; // stop
            }
          break;

        case 5 :   if( DTARGET_print + DELDT > DTARGET0 || EN <= 30. )  {
              DTARGET = DTARGET0;
              flagToStop=false; // stop
            }
          break;
        };
    }

  if(flagToStop)
    goto L5;


  if( optionFlagNorm  &&  iLoop == 5 && fp)
    fprintf(fp,"\n\n  Qmean max = %.2f",QMEANMX);


  if(fp) fclose(fp);

  double qsum=0;
  for(int I=0; I<IMAX; I++)
    {
      if(X[I]<0)                 // v.15.16.7  04/13/21
        {
          if(qsum>0.9) break;
        }
      else
        {
          qsum += X[I];
          Obmen[I]=X[I];
        }
    }

  Obmen[IMAX  ] = QMEAN;
  Obmen[IMAX+1] = dQ;
  Obmen[IMAX+2] = D_EQ;
  Obmen[IMAX+3] = (optionFlagEfinal ? E_init : EOUT);

  return Result;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int FillArrayFromFile(const char *fname, double *ar)
{
  char CFILE_ALL[50];
  char buff[SizeBuf], *t;
  sprintf(CFILE_ALL,"charge_data\\cst%s.data",fname);
  FILE *f9=fopen(CFILE_ALL,"rt");
  if(!f9) return -1;

  for(int ID=0; ID<40; ID++) {
      if(!fgets(buff, SizeBuf,f9)){fclose(f9); return -2;}
      t=buff;
      for(int i=0; i<4; i++) t=GetDoubleFromString(ar[i+ID*4], t);
    }

  fclose(f9);

  return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int GlobalExtern(double *qtab, double *Inform, int fast, double gAF, double gZF, double gAT,
		     double gZT, double gDTARGET, double gEN0, int gQIN, int iOption)
{
  double Obmen[IMAX+4];    //   Qmean, dQ, EquilibThick, Eout
  for(int i=0; i<4; i++) Inform[i]=0;

  for(int i=0; i<IMAX+4; i++) Obmen[i]=0;
  if(gQIN==777) for(int i=0; i<qMin(IMAX,int(gZF));i++) Obmen[i]=qtab[i];

  int result = RunGlobalLocal((double*)Obmen, nullptr, false,
                              fast, nullptr,
                              gAF, gZF, gAT, gZT,
                              gDTARGET, gEN0, gQIN,
                              4, iOption, 2,
                              0, 10, 1,
                              2, 30.,
                              2, 2, 100.);


  if(result >= 0) {
      for(int i=0; i<qMin(IMAX,int(gZF));i++) qtab[i]  = Obmen[i];
      for(int i=0; i<4; i++)                  Inform[i]= Obmen[IMAX+i];
    }

  return result;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
