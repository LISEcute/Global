#include <QtMath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "globallib.h"
#include "gl_data_cst.h"     /* scaled binding energies interpolated ar_fk/l[] */

extern  double ETATH[41], THETA[5];  // increased by 1
extern  double BSA, GA, BETA, BSQ;

// without this double prototyping (mocadi.h and here) type is always wrong
double ffintpol(double ETHH, double THE,
               double TH1,  double TH2,
               double ETH1, double ETH2,
               double F1, double F2, double F3);
double frint1lin(double X, double X1, double X2, double Y1, double Y2);
double frint1log(double X, double X1, double X2, double Y1, double Y2);
double frint1d ( double X, double *XTAB, double *YTAB, int N );
double frint1dl( double X, double *XTAB, double *YTAB, int N );


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void intpoklm(double ETHH, double THE, double Z, double U, const char *CFILE,
              double *SI, double *RTR, double *PP, int *IRC)
{
//
//     interpolates KLM cross sections with F functions
//
      double /*F[4][40]*/ F[5][41], *Fhelp, TH1, TH2, ETH1, ETH2;
      double F1, F2, F3, F4, FTH, BL, CL, C1, C2, C3, D_ETA, DL, ETA; //HWE
      double SLO, HELP1[3], HELP2[3], ETA_T, F_MAX, F_SMALL, AL, EL;
      double F_SMALL_T, FL, G_ETA, G_ETA_T,RL;
      int ID;
// ****************************************************************;
//      double         PI, U0, MC2, S0, ALPH, GA, BSQ, BETA, BSA, D1, D2, GAG, ETA;
//      common /CONST/ PI, U0, MC2, S0, ALPH, GA, BSQ, BETA, BSA, D1, D2, GAG, ETA;
// ****************************************************************;
//      common /ETATH/  ETATH, THETA;
// ****************************************************************;
      int flag=0;
      int i,i1,i2 ;
//
      double FTHE1[5], FTHE2[5]; // increased on 1

     // increased on 1

      double THE_TAB[23] ={0, 1.0, .95, .90, .85, .80, .75, .70, .68, .66, .64, .62,   .60, .58, .56, .54, .52, .50, .48, .46, .44, .42, .40 };
      double C1_TAB[23] ={ 0, 0.2834, 0.3264, 0.3784, 0.4422, 0.5211, 0.6200, 0.7458, 0.8056, 0.8721, 0.9462, 1.0290, 1.1218, 1.2262, 1.3442, 1.4779, 1.6303, 1.8046, 2.0050, 2.2369, 2.5065, 2.8221, 3.1943 };
      double C2_TAB[23] ={ 0, 1.6476, 1.6806, 1.6972, 1.6887, 1.6427, 1.5404, 1.3538, 1.2465, 1.1154, 0.9560, 0.7630, 0.5299, 0.2488,-0.0899,-0.4981,-0.9905,-1.5852,-2.3054,-3.1798,-4.2455,-5.5499,-7.1547 };
      double C3_TAB[23] ={ 0, 8.0519, 7.8541, 7.6611, 7.4946, 7.3903, 7.4064, 7.6383, 7.8236, 8.0813, 8.4285, 8.8861, 9.4803,10.2433,11.2154,12.4471,14.0021,15.9609,18.4262,21.5302,25.4435,30.3892,36.6612 };
      double ATAB[4]   ={ 0,22.222, 16., 4.028};
      double BTAB[4]   ={ 0, 6.333, 8.667, 11.469};
      double CTAB[4]   ={ 0, 4.741, 12.642, 3.330};
// *********************************************************************
//
      *IRC = -1;
      *SI = *RTR = *PP = 0.;

      if( ETHH == 0. || U == 1. ) return;

      for(i1=0; i1<=4; i1++)
            for(i2=0; i2<=40; i2++) F[i1][i2]=0;

      int MODE=-1;

           if( !strncmp(CFILE,"fk", 2)) {  MODE=0; Fhelp=ar_fk0; }
      else if( !strncmp(CFILE,"fl1",3)) {  MODE=1; Fhelp=ar_fl1; }
      else if( !strncmp(CFILE,"fl2",3)) {  MODE=1; Fhelp=ar_fl2; }
      else if( !strncmp(CFILE,"fm1",3)) {  MODE=2; Fhelp=ar_fm1; }
      else if( !strncmp(CFILE,"fm2",3)) {  MODE=2; Fhelp=ar_fm2; }
      else if( !strncmp(CFILE,"fm3",3)) {  MODE=2; Fhelp=ar_fm3; }

      if(MODE < 0) return;



//  sprintf(CFILE_ALL,"charge_data\\cst%s.data",CFILE);
//  f9=mfopen(CFILE_ALL,"rt");

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//     Next line starts K shell interpolation
//L10:
if(MODE==0) {
//====================================
	if( ETHH < 0.01 ) {
                double B0 = 2912.711;
                double B1 = 19.* THE - 43.636;
                double B2 = 102.857 * ( 2.091*pow2(THE)- 10.778*THE + 12.308 );
                FTH = B0 *pow_int(ETHH,4)* ( 1. + B1*ETHH + B2*pow2(ETHH));
                FTH /= THE;
        }
        else if( ETHH >= 0.01 && ETHH <= 86.6 ) {

               for(ID=1; ID<=40; ID++) {

                    for(i=1; i<=4; i++)
                        F[i][ID] = Fhelp[(ID-1)*4 + (i-1)];

                    if( ETATH[ID] > ETHH && ID >= 2 ) break;
               }

               if( THE <  THETA[2])                          flag=0;
               else if( THE >= THETA[2] && THE <= THETA[3] ) flag=1;
               else                                          flag=2;

               TH1=THETA[1+flag];
               TH2=THETA[2+flag];
               ETH1=ETATH[ID-1];
               ETH2=ETATH[ID];
               F1=F[1+flag][ID-1];
               F2=F[1+flag][ID];
               F3=F[2+flag][ID-1];
               F4=F[2+flag][ID];   // at TH2, ETH2   HWE

	       F1 = frint1lin(ETHH,ETH1,ETH2,F1,F2);  // HWE
	       F2 = frint1lin(ETHH,ETH1,ETH2,F3,F4);  // HWE
	       FTH = frint1log(THE,TH1,TH2,F1,F2)/THE;// HWE
	       
               //printf(" THE=%f %f ETH=%f %f F=%f %f %f %f\n",TH1,TH2,ETH1,ETH2,F1,F2,F3,F4); //HWE
	       //printf("  -K %d THE=%f ETHH=%f FTH=%f %f %f \n",ID,THE,ETHH,FTH*THE,F1,F2); // HWE
               //FTH = ffintpol(ETHH,THE,TH1,TH2,ETH1,ETH2,F1,F2,F3); //HWE old
	       //printf("  -K %d THE=%f ETHH=%f FTH=%f \n",ID,THE,ETHH,FTH*THE); // HWE
	       
       }
      else {
             for(i=4; i>=1; i--) FTHE1[i] = Fhelp[39*4 + (4-i)];


                C1 = frint1d(THE,THE_TAB,C1_TAB,22);
                C2 = frint1d(THE,THE_TAB,C2_TAB,22);
                C3 = frint1d(THE,THE_TAB,C3_TAB,22);

                F_MAX = frint1dl(THE,THETA,FTHE1,4);

                ETA_T = 86.6 *pow2(THE);
                F_SMALL_T = F_MAX * ETA_T / THE;
                G_ETA_T = C2 / 4. / ETA_T + C3 / 32. /pow2(ETA_T);

                double log_eta_t =  (ETA_T>0 ? log(ETA_T) : 0);

                D_ETA = F_SMALL_T - C1 * log_eta_t + G_ETA_T;
                ETA = ETHH *pow2(THE);

                G_ETA = C2 / 4. / ETA + C3 / 32. /pow2(ETA);

                double log_eta =  (ETA>0 ? log(ETA) : 0);

                F_SMALL = C1 * log_eta + D_ETA - G_ETA;
                FTH = F_SMALL / ETA;
             }
        }

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//     Next line starts L shell interpolation
else if(MODE==1) {
//=============================
//L40:
    if( ETHH < 0.01 ) {
        for(i=4; i>=1; i--) FTHE1[i] = Fhelp[    (4-i)];
        for(i=4; i>=1; i--) FTHE2[i] = Fhelp[4 + (4-i)];

        HELP2[1] = frint1dl(THE,THETA,FTHE1,4);
        HELP2[2] = frint1dl(THE,THETA,FTHE2,4);
        HELP1[1] = ETATH[1];
        HELP1[2] = ETATH[2];

        FTH  = frint1dl(ETHH,HELP1,HELP2,2);
        FTH /= THE;
  }

  else if( ETHH >= 0.01 && ETHH <= 86.6 ) 
        {
        for(ID=1; ID<=40; ID++) 
               {
                    for(i=1; i<=4; i++)
                        F[i][ID] = Fhelp[(ID-1)*4 + (i-1)];

                    if( ETATH[ID] > ETHH && ID >= 2 ) break;
               }

               if( THE <  THETA[2])                          flag=0;
               else if( THE >= THETA[2] && THE <= THETA[3] ) flag=1;
               else                                          flag=2;

               TH1=THETA[1+flag];
               TH2=THETA[2+flag];
               ETH1=ETATH[ID-1];
               ETH2=ETATH[ID];
               F1=F[1+flag][ID-1]; // at TH1, ETH1
               F2=F[1+flag][ID];   // at TH1, ETH2
               F3=F[2+flag][ID-1]; // at TH2, ETH1
	       F4=F[2+flag][ID];   // at TH2, ETH2   HWE

	       F1 = frint1lin(ETHH,ETH1,ETH2,F1,F2);  // HWE
	       F2 = frint1lin(ETHH,ETH1,ETH2,F3,F4);  // HWE
	       FTH = frint1log(THE,TH1,TH2,F1,F2)/THE;// HWE

         }

  else   {
           for(i=4; i>=1; i--) FTHE1[i] = Fhelp[38*4 + (4-i)];
           for(i=4; i>=1; i--) FTHE2[i] = Fhelp[39*4 + (4-i)];


           HELP2[1] = frint1dl(THE,THETA,FTHE1,4);
           HELP2[2] = frint1dl(THE,THETA,FTHE2,4);
	   HELP1[1] = ETATH[39];
           HELP1[2] = ETATH[40];
           FTH = frint1dl(ETHH,HELP1,HELP2,2);
           FTH /= THE;
         }

   }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//
//
else {   // MODE==2

//     Next line starts M shell interpolation
//L70
      if( ETHH < .04 ) {


             if(!strncmp(CFILE,"fm1",3)) flag=0;
        else if(!strncmp(CFILE,"fm2",3)) flag=1;
        else if(!strncmp(CFILE,"fm3",3)) flag=2;

        RL = flag;
        CL = CTAB[1+flag];
        BL = BTAB[1+flag];
        AL = ATAB[1+flag];

        DL = 1024. * 387420489. *pow(2,2.*RL)*pow(3.,4.*RL)/(5.+RL)/(9.+2.*RL)*CL;
        EL =  324. * (5.+RL) / (6.+RL) * (9.+2.*RL) / (11.+2.*RL) * AL;
        FL =   36. * (5.+RL) / (6.+RL) * (9.+2.*RL) / (10.+2.*RL) * BL;

        FTH = DL *pow(ETHH,4.+RL)* (1.-(EL-FL*THE)*ETHH);
        FTH = FTH / THE;
      }
//-----------------------------------------------------
      else  if( ETHH >= 0.04 && ETHH <= 100. ) {
 
          for(ID=1; ID<=40; ID++) {
                    for(i=4; i>=1; i--)
                        F[i][ID] = Fhelp[(ID-1)*4 + (4-i)];

                    if( ETATH[ID] > ETHH && ID >= 2 ) break;
	  }

          if( THE <  THETA[2])                          flag=0;
          else if( THE >= THETA[2] && THE <= THETA[3] ) flag=1;
          else                                          flag=2;

          TH1=THETA[1+flag];
          TH2=THETA[2+flag];
          ETH1=ETATH[ID-1];
          ETH2=ETATH[ID];
          F1=F[1+flag][ID-1]; // at TH1, ETH1   
          F2=F[1+flag][ID];   // at TH1, ETH2
          F3=F[2+flag][ID-1]; // at TH2, ETH1
	  F4=F[2+flag][ID];   // at TH2, ETH2   HWE
	  
	  //printf(" THE=%f %f ETH=%f %f F=%f %f %f %f\n",TH1,TH2,ETH1,ETH2,F1,F2,F3,F4); //HWE
          F1 = frint1lin(ETHH,ETH1,ETH2,F1,F2);   // HWE
	  F2 = frint1lin(ETHH,ETH1,ETH2,F3,F4);   // HWE
	  FTH = frint1log(THE,TH1,TH2,F1,F2)/THE; // HWE
	  //printf("  -M %d THE=%f ETHH=%f FTH=%f %f %f \n",ID,THE,ETHH,FTH*THE,F1,F2); //HWE
          //FTH = ffintpol(ETHH,THE,TH1,TH2,ETH1,ETH2,F1,F2,F3); // HWE old
	  //printf("  -M %d THE=%f ETHH=%f FTH=%f old \n",ID,THE,ETHH,FTH*THE); //HWE

     }
//-----------------------------------------------------------
     else {
          for(i=4; i>=1; i--) FTHE1[i] = Fhelp[38*4 + (4-i)];
          for(i=4; i>=1; i--) FTHE2[i] = Fhelp[39*4 + (4-i)];

          HELP2[1] = frint1dl(THE,THETA,FTHE1,4);
          HELP2[2] = frint1dl(THE,THETA,FTHE2,4);
	  HELP1[1] = ETATH[39];
          HELP1[2] = ETATH[40];
          FTH = frint1dl(ETHH,HELP1,HELP2,2);
          FTH /= THE;
     };

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//     longitudinal cross.

//      fclose(f9);

      if( FTH <= 0. ) FTH = 0.;

      SLO=7.038E8*FTH/pow_int(Z,4);
//     transverse/longitudinal

      double log_GA = (GA>0 ? log(GA) : 0);
      *RTR=(2.*log_GA-pow2(BETA))/log(1022000.*pow2(BETA)/U);
      if(1022000.*pow2(BETA)/U<1.15) {
        /* quick fix against resonance behavior, correction should always lead to increase HWE */
        *RTR = 0. ;
      }
      *SI=SLO*(1.+ *RTR);
      *IRC = 0;
}



//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double ffintpol(double ETHH, double THE,
               double TH1,  double TH2,
               double ETH1, double ETH2,
               double F1,   double F2, double F3 )
{
//
//     TWO-DIMENSIONAL INTERPOLATION OF F
// This function must be quite wrong, in any case one should use 4 points not only 3, HWE
// **********************************************************************
      if( ETH1 == 0.  || ETH2 == 0. ) return 0;

      double LT1=TH1  > 0 ? log(TH1) : 0;
      double LT2=TH2  > 0 ? log(TH2) : 0;
      double LE1=ETH1 > 0 ? log(ETH1) : 0;
      double LE2=ETH2 > 0 ? log(ETH2) : 0;

      double R12=F1*ETH1/F2/ETH2;
      double B=(R12*LE2-LE1)/(1.-R12)-2.*LT1;
      double L11=LE1+2.*LT1+B;
      double L12=LE1+2.*LT2+B;

      if(  F3*L11/F1/L12 <= 0. ){   
         //printf("  error! expr= %g  f3=%g  L11=%g  F1=%g L12=%g LT1=%g LT2=%g B=%g\n",F3*L11/F1/L12, F3, L11, F1, L12, LT1, LT2, B);
         return 0 ;
      }

      double x = F3*L11/F1/L12;
      double log_x = x > 0 ? log(x) : 0;
      double PP=log_x/(LT1-LT2);

	  TH1 = qMax(1e-50,TH1);
	  THE = qMax(1e-50,THE);

      double A=F1*ETH1*pow(TH1,PP)/L11;

      //     this is F(Benka,Johnson)/THETA
      x = ETHH*THE*THE;
      log_x = x> 0 ? log(x) : 0;

      double FINTPOL=A*(log_x+B)/ETHH/pow(THE,(PP+1.));
      
      return FINTPOL;
}
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// new interpolation functions, old ones had error, H. Weick Jan.2017

/* a simple linear inter (or extra-)polation, HWE */
double frint1lin(double X, double X1, double X2, double Y1, double Y2)
{
  double FRINT ;
  if (fabs(X1-X2)>1e-8)
     FRINT = Y1 +(Y1-Y2)/(X1-X2)*(X-X1);
  else
     FRINT = (Y1+Y2)/2. ;
  
  return FRINT ;
}

/* a linear inter (or extra-)polation on log,log-scale, HWE */
double frint1log(double X, double X1, double X2, double Y1, double Y2)
{
  double FRINT,a,b ;
  if (fabs(X1-X2)>1e-8) {
     X1 = log(X1);
     X2 = log(X2);
     Y1 = log(Y1);
     Y2 = log(Y2);
     b = (Y1-Y2)/(X1-X2);
     a = Y1 - X1 * b;
     FRINT = exp(a) * pow(X,b);     
  }
  else
     FRINT = (Y1+Y2)/2. ;

  return FRINT ;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double frint1d( double X, double *XTAB, double *YTAB, int N )
{
//
//     INTERPOLATION FOR LINEAR FUNCTIONS
//
      int I, ISAVE=1;

// *********************************************************************
      for(I=1; I<=N; I++)
        if( X == XTAB[I] )  return YTAB[I];

           if( X < XTAB[1] ) ISAVE = 1;
      else if( X > XTAB[N] ) ISAVE = N - 1;
      else
            for(I=1; I<=N; I++)
                  if( X < XTAB[I] ) {ISAVE = I - 1; break;}


       double XX = X - XTAB[ISAVE];
       double v = XTAB[ISAVE+1] - XTAB[ISAVE];
       if(fabs(v)<1e-20)v=1e-20;

       double Y = YTAB[ISAVE] + XX * (YTAB[ISAVE+1] - YTAB[ISAVE]) /v;

return Y;
};

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double frint1dl( double X, double *XTAB, double *YTAB, int N )
{
//
//     INTERPOLATION FOR LOGARITHMIC FUNCTIONS
//
      int I, ISAVE=1;
//
// *********************************************************************
//char *str_mess = "\n<E> RINT1DL: NEGATIVE VALUES NOT ALLOWED!\n";

      if( X <= 0.) X = 1.E-50;
      if( X < 0. ) {
         /* fprintf(fp,str_mess);*/  return 0;  };

      for(I=1; I<=N; I++)
        {
           if( XTAB[I] < 0. || YTAB[I] < 0. ) {
               return 0;
           }
           if( XTAB[I] == 0. ) XTAB[I] = 1.E-50;
           if( YTAB[I] == 0. ) YTAB[I] = 1.E-50;
	   if( X == XTAB[I] ) return YTAB[I];
      }


     if( X < XTAB[1] ) ISAVE = 1;
     else if( X > XTAB[N] ) ISAVE = N - 1;
     else
         for(I=1; I<=N; I++)
                if( X < XTAB[I] ) {
                        ISAVE = I - 1;
                        break;
                }

      double l_x    =  (X>0             ? log(X) : 0);
      double l_xis  =  (XTAB[ISAVE]  >0 ? log(XTAB[ISAVE])   : 0);
      double l_xis1 =  (XTAB[ISAVE+1]>0 ? log(XTAB[ISAVE+1]) : 0);
      double XX = l_x  - l_xis;



      double v = l_xis1 - l_xis;

       if(fabs(v)<1e-20)v=1e-20;

      double l_yis  =  (YTAB[ISAVE]  >0 ? log(YTAB[ISAVE])   : 0);
      double l_yis1 =  (YTAB[ISAVE+1]>0 ? log(YTAB[ISAVE+1]) : 0);

      double Y = l_yis + XX  * (l_yis1 - l_yis) / v;
      double RINT1DL = exp(Y);
return RINT1DL;
};

