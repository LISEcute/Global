extern double pow2(double par);
extern double pow_int(double par, int power);



void global_cross(int i_zf, int ne, int i_zm, double EN, double *d_s);
void global_cross2(int i_zf, int i_zm, double EN, double d_s[][4]);

/* gl_ionicro.c (ionization from GLOBAL) */
void ionicro(double ZF, double U1S, double U1SS, double U2S,
             double U2P, double U3S, double U3D, int *IRC);


/* gl_capcross.c (electron captute from GLOBAL) */
void capcross(double ZP, double ZT, double *SC, int *IRC );
void reccork(double Z,  double U, double *SRKBEST, int *IRC);
void reccorlm(double U, double *SRL, double *SRM, int *IRC );


/* gl_bipoco.c (from GLOBAL) */
void bipoco1s(double X, double CX, double ETH, double THE, double Z, double *BP);
void bipoco2s(double X, double CX, double ETH, double THE, double Z, double *BP);
void bipoco2p(double X, double CX, double ETH, double THE, double Z, double *BP);
double fnpol(double X);

/* gl_intpoklm.c (interpolation from GLOBAL) */
void intpoklm(double ETHH, double THE, double Z, double U, const char *CFILE,
              double *SI, double *RTR, double *PP, int *IRC);

double ffintpol(double ETHH, double THE,
               double TH1,  double TH2,
               double ETH1, double ETH2,
               double F1, double F2, double F3);
double frint1d ( double X, double *XTAB, double *YTAB, int N );
double frint1dl( double X, double *XTAB, double *YTAB, int N );
double frint1lin(double X, double X1, double X2, double Y1, double Y2);
double frint1log(double X, double X1, double X2, double Y1, double Y2);


void global_eqdist3(int zp,double ein,int zt,double *f);
void global_eqdist(int zp,double ein,int zt,double *f);
double global_sigma(int zp, double ein, int zt,int ni, int no);
double global_qmean(int i_zf,int i_zt,double e);

/* gl_cross.c  (GLOBAL) */
void global_cross(int i_zf,            /* projectile */
                  int i_qf,            /* number of electrons */
                  int i_zm,            /* target */
                  double d_energy,     /* energy MeV/u */
                  double d_s[]);       /* array of cross sections */

#ifndef PI
#define      PI    3.14159265
#endif

#define      U0    13.6
#define      MC2   511006.
#define      S0    7.038E+08
#define      ALPH  7.2993E-3
#define      IMAX  28

enum CHAR_NAME {
opt_NormIni,
opt_NormEnd,
opt_EqIni,
opt_EqEnd,
opt_Evolution
};
