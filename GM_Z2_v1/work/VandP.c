#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */
extern int access(const char *pathname, int mode);

int nModelParticles=23;
static ModelPrtclsStr ModelPrtcls_[23]=
{
  {"ve","ve~",0, 12, "0","0",1,1,2,0}
, {"vm","vm~",0, 14, "0","0",1,1,2,0}
, {"vt","vt~",0, 16, "0","0",1,1,2,0}
, {"e-","e+",0, 11, "0","0",1,1,2,-3}
, {"m-","m+",0, 13, "MM","0",1,1,2,-3}
, {"tt-","tt+",0, 15, "MTA","0",1,1,2,-3}
, {"u","u~",0, 2, "0","0",1,3,6,2}
, {"c","c~",0, 4, "0","0",1,3,6,2}
, {"t","t~",0, 6, "MT","WT",1,3,6,2}
, {"d","d~",0, 1, "0","0",1,3,6,-1}
, {"s","s~",0, 3, "0","0",1,3,6,-1}
, {"b","b~",0, 5, "MB","0",1,3,6,-1}
, {"A","A",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","WZ",2,1,3,0}
, {"W+","W-",0, 24, "MW","WW",2,1,3,3}
, {"G","G",1, 21, "0","0",2,8,16,0}
, {"h","h",1, 25, "Mh","Wh",0,1,1,0}
, {"~H","~H",1, 252, "MH","WH",0,1,1,0}
, {"H3p","H3p~",0, 253, "MH3","WH3p",0,1,1,3}
, {"~H3z","~H3z",1, 254, "MH3","WH3z",0,1,1,0}
, {"H5pp","H5pp~",0, 255, "MH5","WH5pp",0,1,1,6}
, {"H5p","H5p~",0, 256, "MH5","WH5p",0,1,1,3}
, {"~H5z","~H5z",1, 257, "MH5","WH5z",0,1,1,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=22;
int nModelFunc=14;
static int nCurrentVars=21;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[36]={
 "cabi","lam2","lam3","lam4","lam5","mu3sq","aEWM1","Gf","aS","tanth"
,"sa","ymb","ymt","ymtau","MM","MTA","MT","MB","MZ","Mh"
,"E","Pi","ca","CKM1x1","CKM1x2","CKM2x1","CKM2x2","sh","aEW","v"
,"vphi","ch","MH","MH3","MH5","MW"};
char**varNames=varNames_;
static REAL varValues_[36]={
   2.277360E-01,  1.000000E-01,  0.000000E+00,  1.000000E-01,  1.000000E-01,  1.000000E-01,  1.279000E+02,  1.166370E-05,  1.184000E-01,  0.000000E+00
,  0.000000E+00,  4.700000E+00,  1.720000E+02,  1.777000E+00,  1.056600E-01,  1.777000E+00,  1.720000E+02,  4.700000E+00,  9.118760E+01,  1.250000E+02
,  2.718282E+00,  3.141593E+00};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)     return 0;
   }
  cErr=1;
   nCurrentVars=22;
   V[22]=Pow(1-Pow(V[10],2),0.5);
   if(!isfinite(V[22]) || FError) return 22;
   nCurrentVars=23;
   V[23]=Cos(V[0]);
   if(!isfinite(V[23]) || FError) return 23;
   nCurrentVars=24;
   V[24]=Sin(V[0]);
   if(!isfinite(V[24]) || FError) return 24;
   nCurrentVars=25;
   V[25]=-Sin(V[0]);
   if(!isfinite(V[25]) || FError) return 25;
   nCurrentVars=26;
   V[26]=Cos(V[0]);
   if(!isfinite(V[26]) || FError) return 26;
   nCurrentVars=27;
   V[27]=V[9]*Pow(1+Pow(V[9],2),-0.5);
   if(!isfinite(V[27]) || FError) return 27;
   nCurrentVars=28;
   V[28]=Pow(V[6],-1);
   if(!isfinite(V[28]) || FError) return 28;
   nCurrentVars=29;
   V[29]=Pow(2,-0.25)*Pow(V[7],-0.5);
   if(!isfinite(V[29]) || FError) return 29;
   nCurrentVars=30;
   V[30]=Pow(2,0.25)*Pow(V[7],0.5);
   if(!isfinite(V[30]) || FError) return 30;
   nCurrentVars=31;
   V[31]=Pow(1-Pow(V[27],2),0.5);
   if(!isfinite(V[31]) || FError) return 31;
   nCurrentVars=32;
   V[32]=Pow(V[5]+2*V[1]*Pow(V[30],2)-V[4]*Pow(V[30],2),0.5);
   if(!isfinite(V[32]) || FError) return 32;
   nCurrentVars=33;
   V[33]=Pow(2,-0.5)*Pow(2*V[5]+4*V[1]*Pow(V[30],2)-V[4]*Pow(V[30],2),0.5);
   if(!isfinite(V[33]) || FError) return 33;
   nCurrentVars=34;
   V[34]=Pow(2,-0.5)*Pow(2*V[5]+4*V[1]*Pow(V[30],2)+V[4]*Pow(V[30],2),0.5);
   if(!isfinite(V[34]) || FError) return 34;
   nCurrentVars=35;
   V[35]=Pow(Pow(V[18],2)/(2.)+Pow(Pow(V[18],4)/(4.)-V[28]*V[21]*Pow(2,-0.5)*Pow(V[7],-1)*Pow(V[18],2),0.5),0.5);
   if(!isfinite(V[35]) || FError) return 35;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
