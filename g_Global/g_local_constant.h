#ifndef LOCAL_CONSTANT_H
#define LOCAL_CONSTANT_H
//#define      PI    3.14159
#define      U0    13.6
#define      MC2   511006.
#define      S0    7.038E+08
#define      ALPH  7.2993E-3
#define      IMAX  28

#define pow(x,y)   exp((y)*log(x))
#define max(x,y)   (x > y ? x : y)
#define min(x,y)   (x < y ? x : y)

//extern double pow2(double par);
//extern double pow_int(double par, int power);


enum CHAR_NAME {
opt_NormIni,
opt_NormEnd,
opt_EqIni,
opt_EqEnd,
opt_Evolution
};
#endif // LOCAL_CONSTANT_H
