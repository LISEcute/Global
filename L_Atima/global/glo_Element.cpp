#include <cstring>
#include "string_utils.h"

void ELEMENT(double &Z, double &A, char *CZ, int IOPT, int &IRC );
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ELEMENT(double &Z, double &A, char *CZ, int IOPT, int &IRC ) {
//
//     RETURNS ELEMENT SYMBOL FROM NUCLEAR CHARGE
//
//     IOPT 1: CHARGE NUMBER TO ELEMENT SYMBOL CONVERSION
//     IOPT 2: ELEMENT SYMBOL TO CHARGE NUMBER CONVERSION
//     IOPT 3: CHARGE NUMBER TO ELEMENT SYMBOL CONVERSION
//     IOPT 4: ELEMENT SYMBOL TO CHARGE NUMBER CONVERSION
//     IOPT 3  ET 4: INTEGER CUTTING FOR BEAMS
//     IOPT 5: RETURNS Z AND CZ FOR A GIVEN MASS NUMBER A
//     IOPT 6: Z TO SYMBOL CONVERSION
//
//
// **********************************************************************
      const char *CZDATA[98] =
      {"n",  "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ",
       "NE", "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ",
       "CA", "SC", "TI", "V ", "CR", "MN", "FE", "CO", "NI", "CU",
       "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y ",
       "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN",
       "SN", "SB", "TE", "J ", "XE", "CS", "BA", "LA", "CE", "PR",
       "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM",
       "YB", "LU", "HF", "TA", "W ", "RE", "OS", "IR", "PT", "AU",
       "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC",
       "TH", "PA", "U ", "NP", "PU", "AM", "CM", "MY"};
//
      const double MASS[98]={
       1.0,   1.0,    4.0,  6.94,   9.01, 10.81, 12.01, 14.01, 16.00, 19.00,
       20.18,  23.00, 24.31, 26.98, 28.09, 30.97, 32.06, 35.45, 39.95,  39.10,
       40.08, 44.96, 47.88, 50.94, 52.00, 54.94, 55.84, 58.93,  58.69,  63.55,
       65.40, 69.72, 72.59, 74.92, 78.96, 79.90, 83.80,  85.47, 87.62, 88.91,
       91.22, 92.91, 95.93, 98.  , 101.07, 102.91, 106.42, 107.87, 112.41,
       114.82, 118.71, 121.76, 127.60, 126.90, 131.29, 132.91, 137.33,
       138.91, 140.12, 140.91, 144.24, 147., 150.36, 151.97, 157.25, 158.92,
       162.50, 164.93, 167.26, 168.93, 173.03, 174.97, 178.49, 180.95, 183.85,
       186.21, 190.20, 192.22, 195.08, 196.97, 200.60, 204.38, 207.20, 208.98,
       203., 209., 211., 212., 213., 222., 232., 231., 238., 237., 244., 243.,
       247., 13.72};
//
// **********************************************************************
//
      int I;
      if(strlen(CZ)>2) CZ[2]='\x0';
      struprLoc(CZ);  // v.16.18.1
//      strupr(CZ);
      IRC = -1;
//
      switch(IOPT)
        {
        case 1 :
        case 3 :
                if( Z == 6.6 ) { strcpy(CZ,CZDATA[97]);  IRC=0; }
                else
                     for(I=1; I<=96; I++)
                          if( int(Z+0.5) == I ) {
                        strcpy(CZ,CZDATA[I]);
                                IRC=0;
                                break;
                                }
                break;

      case 2 :
      case 4 :
                for(I=1; I<=97; I++)
                        if( strcmp(CZ,CZDATA[I])==0) {
                                    if( strcpy(CZ,"MY")==0) Z = 6.6;
                                    else                    Z = double(I);
                                    IRC = 0;
                                    break;
                                    }
                break;

        case 5 :
                for(I=1; I<=97; I++)
                          if( A <= MASS[I] ) {
                                Z  = double(I);
                                strcpy(CZ,CZDATA[I]);
                                IRC=0;
                                break;
                                }
                break;

      case 6 :  int Zi = int(Z+0.5);  IRC=0;
                if(Zi>=1 && Zi <=97) strcpy(CZ,CZDATA[Zi]);
                else IRC=-1;
                break;
      };


   if(IRC==0) {
      if( Z == 6.6 ) A = MASS[97];
      else           A = MASS[int(Z+0.5)];

      if( IOPT == 3 || IOPT == 4 ) A = int(A+0.5);

      strlwrLoc(&CZ[1]);
      }

}
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//

