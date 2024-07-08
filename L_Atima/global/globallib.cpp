#include "globallib.h"
#include <stdio.h>
#include "math.h"
#include <string.h>


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/**
 * This calculates equilibrium distribution for Zp projectile
 * at the moment distribution only up to 2 electrons are calculated
 */
void global_eqdist(int zp,double ein,int zt,double *f)
{
  double cs[29][4];
  //double ics[29],ccs[29];
  double r[29];
  double psum;
  double d;

  memset(cs,0,sizeof(cs));
  memset(r,0,sizeof(r));
  int max = 28;

#if(0)
  global_cross(zp,0,zt,ein,cs[0]);
  for(int i=1;i<29;i++)
    {
      global_cross(zp,i,zt,ein,cs[i]);
      r[i-1] = cs[i][1]/cs[i-1][2];

      if(r[i-1]>200.){   // threshold to stop calculation higher charge states
          max = i;
          break;
        }
    }
#else
  global_cross2(zp,zt,ein,cs);
  for(int i=1;i<29;i++){
      //r[i-1] = ics[i]/ccs[i-1];

      if(cs[i][1]<=0)   // Oleg v.16.11.6
        {
          max = i;
          r[i-1]=0;
          break;
        }


      r[i-1] = cs[i][1]/cs[i-1][2];

      if(r[i-1]>200.){   // threshold to stop calculation higher charge states
          max = i;
          break;
        }

    }
#endif


  d = 1;
  for(int i=0;i<max;i++){
      psum=1;
      for(int j=i;j<max;j++){
          psum *= r[j];
        }
      d+=psum;
    };

  f[max] = 1/d;
  for(int i=max-1;i>=0;i--){
      f[i] = r[i]*f[i+1];
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/**
 * This calculates equilibrium distribution for Zp projectile
 * distribution only up to 2 electrons are calculated
 */
void global_eqdist3(int zp,double ein,int zt,double *f){
  double cs[3][4];
  double r0,r1,d;

  for(int i=0;i<3;i++){
      global_cross(zp,i,zt,ein,cs[i]);
    }

  if(cs[1][2]<=0 || cs[0][2]<=0 || 0 || (cs[2][1] + cs[1][1]) <=0 )   //v.16.11.6
    {
      r0 = cs[1][1]/cs[0][2];
      r1 = cs[2][1]/cs[1][2];
    }
  else {r0=0; r1=0;}

  d = 1+r1+r0*r1;
  f[0] = r1*r0/d;
  f[1] = r1/d;
  f[2] = 1/d;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double global_qmean(int i_zf, int i_zt, double e)
{
  int i,n ;
  double sum0,sum1/*,sum2*/ ;
  double f[29] ;
  memset(f,0,sizeof(f));
  n=28 ; // max. for GLOBAL

  if((e>=30)&&(e<=2000)&&(i_zf>=29))
    {
      //csd_equilib(e, i_zf, i_zt, n, f);
      global_eqdist(i_zf, e, i_zt, f);

      sum0 = 0. ;
      sum1 = 0. ;

      for(i=0; i<=n; i++)
        {
          sum0 = sum0 + f[i];
          sum1 = sum1 + f[i]*(i_zf-i);
        }

      double div=sum0;
      if(div<=0)
        div=1;
      return  sum1/div;
    }
  else {
      return -1;
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

/**
 * This returns cross section for zp projectile changing from ni to no electrons
 * \param ni incident number of electrons
 * \param no outgoing number of electrons
 * \return cross section ni->no
 */
double global_sigma(int zp, double ein, int zt,int ni, int no)
{
  double cs1[4];
  int dif = no-ni;
  
  if(ni==no)return 0;
  if(fabs(dif)>2)return 0;
  global_cross(zp,ni,zt,ein,cs1);
  //printf("cs = %lf %lf %lf %lf\n",cs1[0],cs1[1],cs1[2],cs1[3]);
  //printf("index = %d\n",(dif>0)?1+dif:2+dif);
  return cs1[(dif>0)?1+dif:2+dif];
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

