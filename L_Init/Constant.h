#if !defined(constant_h)
    #define  constant_h


# define Number_CrossSections 10000      //TO
# define Number_tab 		60000      //TO        modified  12.0.3    11/29/19   it was 40000
# define Number_tab_MC 		30000      //TO     // it was 200 in 9.1.8
# define Number_tab_SR	      20000      //TO        modified  12.0.3    11/29/19   it was 15000
# define Number_tab_iso		5000       // isomer
# define Number_EPAX		6	     // 04/27/2012 - 5    // 04/28/2017 - 6
# define Number_EPAX_SR		5	     // 08/28/2012 -3 ;  05/02/2017 - 4 ; 03/07/2024 - 5
# define Number_Fusion		1
# define NumberTabELOSS		100
# define NumberTabCHARGE      100
# define Number_ELOSS		5      // 11/01/20 + Atima 1.4
# define Number_ELOSS_PLOT	(Number_ELOSS+2)
//# define Number_RangePoints  217
# define Number_RangePoints  268
# define Number_StragEne	3
# define Number_StragAng	4
# define Number_EqThick 	2   //SRIM does wwork very bad
# define Number_LISE_web	4
# define Number_LISE_version_web 4
# define Number_LISE_ftp	1
# define Number_time_new_versions	8     // Ntnv
# define Number_reactions	   9    // 3 - Abrasion-Fissuion
# define Number_SecondaryTargets 3    // 3 - Secondary targets
//# define Number_TI 17    // Time interval
//# define Number_TI 20    // Time interval                        01/02/2013
//# define Number_TI 21    // Time interval                        03/03/2013
# define Number_TI 22    // Time interval                        03/06/2013   + "ys"

# define Number_T12_alpha	   5
# define Number_T12_spontaneousNPA  4 // NPA
# define Number_T12_spontaneous     6 // 4+qMax + Zagrebaev
# define Number_CatcherMaterials    6
# define NumberConvolutionPlotNames 10
#define Number_EE_models 	5     // excitation energy models

//# define Number_Button_reactions  Number_reactions
# define Number_AxisElementModes 9       // 5 regular, 3sum+NdZ
# define Number_ProfilesPhys	3
# define Number_Profiles	4
# define Number_MassFormula   3
# define Number_MassDataBase  2
# define Number_FisBar		7
# define Number_T12		4
# define Number_T12_compile	4  //04/03/2013 changed from 3 to 4
# define Number_SecR_filtres  4
# define Number_RayFields	15
# define Number_RayLocations	15 // 06/15/2021 changed from 12
# define Nfingers 5

# define Number_BeamShape1	4
# define Number_BeamShape2	4

# define Number_Edipo		3
#define BadValue 		    -777

#define strSize1000           1020
#define strSize4000           4080

#define FlagLogTrue    true
#define FlagLogFalse   false
//#define NumberBLOCKs 		13    // v76
//#define NumberBLOCKs 		14    // v83
//#define NumberBLOCKs 		15    // v83.141
//#define NumberBLOCKs 		16    // v9.4.21
//#define NumberBLOCKs 		17    // v9.6.113
#define NumberBLOCKs 		18    // v9.10.78   Fit block
#define NumberDIfields		11
#define NumberGanilIsomers    500
#define NumberZion	 	110
#define NumberZionArray  	(NumberZion+1)


# define Number_Charges		6      ///  charge state methods

# define Number_Proton		200
# define Number_Neutron		500
# define Number_Nucleon		700
# define Number_NP		7		//16,32,64,128,256,512,1024
# define MaxZtarget		103
# define MaxZtargetHubert     92
#define  NumberMethodSeparation 11   // 8 method + fission + breakup + gamma&residual

//------------ Monte Carlo ---------------------  start
//#define  MC_RB 13                    // MC transmission parameters
//#define  MC_RB 17                    // MC transmission parameters     08/14/07


/*  06/18/2013    v.9.6.38
#define  MC_RB 		18                    // MC transmission parameters     03/10/10  + Envelope
#define  MC_RB_AZQ 	9                     // MC AZQ transmission parameters     04/09/13   Mass,Ad,Z,N, q, ZmQ,AoQ, A-2Q,A-3Q
#define  MC_RB_PARAM 	27  //(MC_RB+MC_RB_AZQ)    // MC all transmission parameters     04/09/13   Mass,Ad,Z,N,q, ZmQ,AoQ, A-2Q,A-3Q
*/

//#define  MC_RB 		17                    // MC transmission parameters      - Velocity
#define  MC_RB 		18                    // MC transmission parameters      - Velocity  + Erho 07/19/2013
#define  MC_RB_real 	20                    // MC_RB + 2 groups (VEL and AZQ)
#define  MC_RB_AZQ 	9                     // MC AZQ transmission parameters
//#define  MC_RB_VEL 	5                     // MC AZQ transmission parameters
#define  MC_RB_VEL 	7                     // MC AZQ transmission parameters    + Beta, Gamma  11/21/2020
#define  MC_RB_PARAM 	34  //(MC_RB + MC_RB_VEL + MC_RB_AZQ)    // MC all transmission parameters + Beta, Gamma  11/21/2020


//#define  MC_RB_files 17              // MC transmission parameters for ray files    04/07/08     no velocity?
//#define  MC_RB_files 20              // MC transmission parameters for ray files     02/28/2010
//#define  MC_RB_files 27              // MC transmission parameters for ray files     04/09/2013 Ad, N,ZmQ,AoQ, A-2Q,A-3Q   07/19/2013  + Erho
//#define  MC_RB_files 28              // MC transmission parameters for ray files     04/09/2013 Ad, N,ZmQ,AoQ, A-2Q,A-3Q   07/19/2013  + Erho   + 07/22/2013  reaction
//#define  MC_RB_files 35              // MC transmission parameters for ray files     04/09/2013 Ad, N,ZmQ,AoQ, A-2Q,A-3Q   07/19/2013  + Erho   + 07/22/2013  reaction   + 5 velocities + 2 Veloc 11/21/20
#define  MC_RB_files 37              // MC transmission parameters for ray files     04/09/2013 Ad, N,ZmQ,AoQ, A-2Q,A-3Q   07/19/2013  + Erho   + 07/22/2013  reaction   + 5 velocities + 2 Veloc 11/21/20 + Zb[m] Zend
#define  MC_plot_from_files 12       // MC transmission parameters to plot in File dialog
#define  Number_MCgates	    4
#define  MC_input        5     // 04/11/2013    09/19/2022
//------------ Monte Carlo ---------------------  stop

#define SR_NpickupShift 2		// for fragmentation -- neutron pickup secondary reactions
#define SR_PpickupShift 1           // for fragmentation -- proton pickup secondary reactions


#define  max_Matter_into_Compound 5
#define  N_DecayString           19   // 12/5/2011   + proton & b+     &&   proton & a     // 12/5/2011    b- & a   ,  12/6/2011  SF & b+, 17: SF & b-    ,      08/02/16   Beta- n  ,    08/09/17 Beta+ Beta+P
#define  N_PartialDecay          10		// + Beta+P
#define  LengthNameCompound      76  // modified v.16.3.34
#define  LengthNameCompoundPlusMass 100
//#define  LengthNameCompound (max_Matter_into_Compound*12+1)  // modified v.14.20.2
//#define  LengthNameCompoundPlusMass (LengthNameCompound+17)

#define  N_sigma_corr		   2
#define  N_sigma_distr		   4
#define  N_velocity_distr	   6
#define  N_mode_separation	   3
#define  NumberDriftModes	   5    // 01/14/2014  - quad - sext + multi + finger

// secondary reactions
//#define  n_draw_SecReact	   5     // and +1 for thickness target
#define  n_draw_SecReact	   6     // and +1 for thickness target + 1 sum of lost  04/20/2024
//#define  coef_decrease_SecReact  4
#define  Threshold_Alpha 	   6.0

#define titreSize        strSize4000
#define filenameCFG_size strSize1000
//# define StrLenFileName		50


//   PHYS

#define  ln_2			0.693147180559945309417
# define PI				3.14159265358979323846
# define PI2			6.28318530717959
# define PI_L			3.1415926535897932384626433832795L
# define vc_mm 			299.792458    //TO  mm/ns
# define vc_cm 			29.9792458    //TO  cm/ns
# define vc_m 			0.299792458   //TO  m/ns
# define vc_m_s 			2.99792458e+8 //TO  m/s
# define vc_cm_s 			2.99792458e+10 //TO  cm/s
# define vc_fm_s 			2.99792458e+23  //TO  fm/s

#define _vc              vc_cm_s       // light velocity in cm/sec
#define _sec             1.             //system value of time sec=1
#define _ns              (1.e-9*_sec)
#define _vcn             (_vc*_ns)           // light velocity in cm/ns


#define RAD 	  		0.017453292519943		// (PI/180.)
#define RAD_L	 	  	0.0174532925199432957692369076848861L
#define GRAD		  	57.2957795130823	// (180./PI)
#define GRAD_L			57.2957795130823208767981548141052L

#define RADM 	  		17.453292519943		// (PI/180.)
#define RAD_LM	 	  	17.4532925199432957692369076848861L
#define GRADM		  	0.0572957795130823	// (180./PI)
#define GRAD_LM			0.0572957795130823208767981548141052L

//-----------------------------------
/*#ifndef  __MATH_H
#define  __MATH_H
// Constants rounded for 21 decimals.
#define M_E         2.71828182845904523536
#define M_LOG2E     1.44269504088896340736
#define M_LOG10E    0.434294481903251827651
#define M_LN2       0.693147180559945309417
#define M_LN10      2.30258509299404568402
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923
#define M_PI_4      0.785398163397448309616
#define M_1_PI      0.318309886183790671538
#define M_2_PI      0.636619772367581343076
#define M_1_SQRTPI  0.564189583547756286948
#define M_2_SQRTPI  1.12837916709551257390
#define M_SQRT2     1.41421356237309504880
#define M_SQRT_2    0.707106781186547524401
#endif
*/
#ifndef  M_PI
#define M_PI        3.14159265358979323846
#endif
//-----------------------------------

#define amu_mg           (1.660538921e-21)   // mg
//#define _joule           1.             //system value of Energy joule=1
#define _eV              1.6021766208e-19   //_joule)
#define _e_charge 	 _eV			 // Coulomb
#define _MeV             (1.e6*_eV)
#define _cm              1.             //system value of length cm=1
#define _fermi           (1.e-13*_cm)
//#define _h_bar           6.58211889e-22   // MeV/s
#define _hbarc            6.58211889e-22   // MeV/s           time -> hbabc / eV =   	6.582119?0-16 s
#define _hbarc_MeV_fermi  197.327053     //  hbarc in MeV*fm

//#define _hc              (_hc_MeV_fermi*_MeV*_fermi)
//#define _e2              _hc
#define FineStructConst  7.297352521e-03       // fine-structure constant      (e/Pc)^2   the square of the ration of the elementary charge to the Planck charge
// or  alpha(FineStructConst) = e^2/ _hbarc
#define e2_MeV_fermi    (_hbarc_MeV_fermi*FineStructConst)


//from MathCad
//    hrc == 197.327 MeV*fermi
//    e2 = hrc / 137               from rutherford for example
//


#define barn            (1.e-24*_cm*_cm)
#define mbarn           (1.e-3*barn)
//#define U               amu_MeV
#define Brho_const	3.107129578
#define Erho_const      931494013.5   // Brho_const*v_m_s   06/02/2024  v.17.5.16
//#define Erho_const      9.32764e8   // Joule/coul  wrong


//from the code "charge"
#define BohrRadius_cm      0.529177249e-08       // Bohr radius in cm;
#define Rydberg_eV         13.605698100000D00    // Rydberg in eV;
// for Schiwietz
#define Bohr_Velocity    2.19e6  // m/s

#define amu_MeV		931.494061  //TO  MeV      12/05/2011 from Wikipedia

#define Mass_p          938.272046          // MeV 12/05/2011 from Wikipedia           nucleon
#define Mass_n          939.565378          // MeV 12/05/2011 from Wikipedia           nucleon

#define Mass_e          0.510998910          // MeV   10/14/2013 from Wikipedia

#define ME_1H            7.2889705          // MeV       atom
#define ME_1n            8.0713171	        // MeV       atom
#define ME_2H            13.1357216	        // MeV       atom
#define ME_3H            14.9498060	        // MeV       atom
#define ME_3He          14.9312148	        // MeV       atom
#define ME_4He           2.42491565	        // MeV       atom

#define DontCalculateCSdegrader   -7711
//-----------------------------------------------------------
#define  o1_div_32767	(3.05185e-5)

#define  local_erf_value  4.1
#define  max_erf_value	  19.9

#define calc_index(Z,N)          (Z*1000+N)
#define calc_indexG(Z,N,G)       ((Z*1000+N)*1000 + int(G))
#define Z_from_index(index)   (index/1000)
#define N_from_index(index)   (index%1000)

#define exp_ln_limit  1e-43 //   exp(-99) =  1.011e-43
#define min_gauss  1e-6

#define DefaultExEnerAA 14.
#define DefaultMassOption 201
#endif
