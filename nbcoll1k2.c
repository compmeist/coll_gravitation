

/*   
    Constrained N-Body simulation to test collateral gravity hypothesis
      (first order test) 
      Nathan E. Frick
*/



/* 

                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "[]"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright 2023 Nathan E. Frick

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/



#include <stdio.h>
#include <math.h>
#include <malloc.h>


#define sqr(x) ((x)*(x))
#define PI 3.14159265358979323846 

/* prototype function definitions */

void initRand(int *idummy );
float randPM (float maxOut,int *idummy); /* plus maxOut to minus maxout */
float randP (float maxOut,int *idummy);
float gaussianDeviate (float mean,float variance, int *idummy);

/* data structures */

typedef struct bdType {
   int id;   /* body number */
   double c1,c2,c3;  /* center point (position of body) */
   double phi,theta,cm;  /* direction (from origin) associated with this position */
   double r,rplus;   /* radius of spherical body */ 
   double m,gx,gy,gz;  /* mass and gravitational accel */
   double cgx,cgy,cgz;  /* collateral gravity */
   double vx,vy,vz;  /* velocity of body */
   double ampF,ampDist;   /* amplification factor */
   double r2,c2_r2,discrim,t,c_dot_d,discrimP;  /* precalc variables */
   double c2_rplus2;
} bdType;

/* computational pages: designed to allow for millions of elements, if needed */

#define PARY_MAX 65520

typedef struct cmpPageNodeType {  /* break the large data into a list of pages (arrays) */
  int status;
  int n;
  bdType ary[PARY_MAX];   /* instead of a monolithic array, take this approach */
  struct cmpPageNodeType *nextnode;
}  cmpNode;

/* define a type cmpNode for convenience */
 
cmpNode cpNode0,*curNode;   /* cpNode0 is the the first page (base node) */
int largestID,obsID;
double collFraction,initVelocity;

char confine2Plane[128];

int allocMem4Arrays(int nbTotal)
/* initialization of memory - returns zero for normal alloc */
/* this modular approach does not require a giant block of memory from the system,
      and lends itself to multithreaded implementation, if required */
{  cmpNode *curNode; int i;
  
   curNode = &cpNode0; /* start at the base node */
   curNode->n = 0;
   for (i=0;i < nbTotal;i++)
   { if (curNode->n == PARY_MAX)
     { /* get next node */
       curNode->nextnode = (cmpNode *)malloc(sizeof(cmpNode));
       if (curNode->nextnode == NULL) return 1;
       curNode = curNode->nextnode;
       curNode->n = 0;
     }
     curNode->ary[curNode->n].id= i;   /* assign an enumeration id */
     curNode->status = 0;
     curNode->n++;
    
   }
   return 0;
}



void CheckForCollionsMerge(void)

/* here is a simple elastic merge of bodies into one */

{ cmpNode *curNode,*jcurNode; int i,j; int curId;
  double Fpi3;  double dx,dy,dz,dist2,r1,r2,totM,totV,averageDensity;
  double totalMomentumX,totalMomentumY,totalMomentumZ;
  double newV1,newM1,newR1,newP1x,newP1y,newP1z;
  Fpi3 = 4 * PI / 3.0;
   curNode = &cpNode0; /*  base node */
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)  if (curNode->ary[i].id != largestID) /* don't allow largest body to merge */
      {  curId = curNode->ary[i].id;
       jcurNode = &cpNode0;
     while (jcurNode != NULL) {
     for (j=0;j<jcurNode->n;j++) if (jcurNode->ary[j].id != curId)
         { dx = jcurNode->ary[j].c1 - curNode->ary[i].c1;
       dy = jcurNode->ary[j].c2 - curNode->ary[i].c2;
       dz = jcurNode->ary[j].c3 - curNode->ary[i].c3;
       dist2 = dx * dx + dy * dy + dz * dz;
       if (dist2 < sqr(jcurNode->ary[j].r +  curNode->ary[i].r) )
       { /* collision */
         r1 =  curNode->ary[i].r;
       r2 =  jcurNode->ary[j].r;
         totM = jcurNode->ary[j].m + curNode->ary[i].m;  /* get various totals */
       totV = Fpi3 * (r1*r1*r1 + r2*r2*r2);
       averageDensity =   totM / totV;
       totalMomentumX =  jcurNode->ary[j].m * jcurNode->ary[j].vx +   curNode->ary[i].m *  curNode->ary[i].vx;
       totalMomentumY =  jcurNode->ary[j].m * jcurNode->ary[j].vy +   curNode->ary[i].m *  curNode->ary[i].vy;
             totalMomentumZ =  jcurNode->ary[j].m * jcurNode->ary[j].vz +   curNode->ary[i].m *  curNode->ary[i].vz;
       newV1 = totV;     /* make some assumptions here regarding volume */
       newM1 = totM;
       newR1 = pow((newV1 / Fpi3),(1/3.0));
       newP1x = totalMomentumX; 
       newP1y = totalMomentumY;
       newP1z = totalMomentumZ;
       if (r1 > r2)  /* let the previously larger body remain so */
       {  curNode->ary[i].r = newR1;
          curNode->ary[i].m = newM1;
        curNode->ary[i].vx = newP1x / newM1;
        curNode->ary[i].vy = newP1y / newM1;
        curNode->ary[i].vz = newP1z / newM1;
        jcurNode->ary[j].r = 0.0;
          jcurNode->ary[j].m = 0.0;
        
       }
       else
       {  curNode->ary[i].r = 0.0;
          curNode->ary[i].m = 0.0;
        jcurNode->ary[j].r = newR1;
          jcurNode->ary[j].m = newM1;
        jcurNode->ary[j].vx = newP1x / newM1;
        jcurNode->ary[j].vy = newP1y / newM1;
        jcurNode->ary[j].vz = newP1z / newM1;
       }
       
         }
     }
             jcurNode = jcurNode->nextnode;
    }
        }
      curNode = curNode->nextnode;
   }   
  
}

void DumpZeroBodies( cmpNode *node0)
{ /* one step in ray-sphere intersection calc */
  cmpNode *curNode; int i,j;
   curNode = node0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
      { if (curNode->ary[i].m == 0.0)
        {  for (j=i+1;j<curNode->n;j++)  curNode->ary[j-1] = curNode->ary[j];  /* shift remaining elements over */
           curNode->n = curNode->n - 1;
        }
      }     
      curNode = curNode->nextnode;
   }
}

void DumpFarBodies( cmpNode *node0,double ex)
{ /* one step in ray-sphere intersection calc */
  cmpNode *curNode; int i,j;
   curNode = node0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
      { if (  (fabs(curNode->ary[i].c1) > ex) ||
              (fabs(curNode->ary[i].c2) > ex) ||
              (fabs(curNode->ary[i].c3) > ex) )
        {  for (j=i+1;j<curNode->n;j++)  curNode->ary[j-1] = curNode->ary[j];  /* shift remaining elements over */
           curNode->n = curNode->n - 1;
        }
      }     
      curNode = curNode->nextnode;
   }
}


void buildBodiesN(double initExtent, double avgDensity, double sigDensity, 
         double smallestRadius, double incRadius,double rFacBump,
         int iRandSkip)

/* seed a normal distribution to the mass density, and a lognormal distribution to the radius */

{  double m,r,rnd,rndRadius,largest,rFac;
   int j,iTotal,iTt,idk;
   cmpNode *curNode; int i; int idmy;
   initRand(&idmy);  /* seed random number gen */
   /* first, seed radii based on a sort of lognormal */
   iTotal = 0;
   curNode = &cpNode0; /*  base node */
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
      { iTotal++;
        curNode->ary[i].r = smallestRadius; 
      }
      curNode = curNode->nextnode;
   }
   rFac = rFacBump;
   iTt = iTotal /  2.0;
   while (iTt > 0)
   { j = 0;
     while (j < iTt)
     { idk = randP(iTotal,&idmy);
       curNode = &cpNode0; /*  base node */
       while (curNode != NULL)
       { for (i=0;i<curNode->n;i++)
         if (idk == curNode->ary[i].id)
        { curNode->ary[i].r += rFac * smallestRadius; 
        }
         curNode = curNode->nextnode;
       } 
       j++;
     }
     iTt /= 2.0;
     rFac *= rFacBump;    /* bump factor for next round of increments */
   }
   curNode = &cpNode0;  /* attach a gaussian noise */
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) curNode->ary[i].r *= (1 + gaussianDeviate(0.0,0.16,&idmy));     
      curNode = curNode->nextnode;
   }

   /* now, seed the masses */
   curNode = &cpNode0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)curNode->ary[i].m = gaussianDeviate(avgDensity,sigDensity,&idmy);    
      curNode = curNode->nextnode;
   }
   curNode = &cpNode0; /* check to make sure there are no negative values */
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) if (curNode->ary[i].m < 0.0) curNode->ary[i].m = avgDensity;      
      curNode = curNode->nextnode;
   }
   
   /* now, multiply the volume and density to get total mass per body */
   curNode = &cpNode0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
      { r = curNode->ary[i].r;
        curNode->ary[i].m = curNode->ary[i].m * (4 * PI / 3) * r * r * r;    
      }
      curNode = curNode->nextnode;
   }
   /* determine which elements is very largest */
   curNode = &cpNode0;
   largest= 0.0;  largestID = 0;
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
      if (curNode->ary[i].m > largest)  {
         largestID = curNode->ary[i].id;
         largest = curNode->ary[i].m;
      }
      curNode = curNode->nextnode;
   }
   /* assign random position within a given space */
   curNode = &cpNode0;  initRand(&idmy);
   for (i=0;i<=iRandSkip;i++) randPM(initExtent,&idmy);  /* seed random num */
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) curNode->ary[i].c1 = randPM(initExtent,&idmy);
      curNode = curNode->nextnode;
   }
   curNode = &cpNode0;  
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) curNode->ary[i].c2 = randPM(initExtent,&idmy);  /* y coordinate */
      curNode = curNode->nextnode;
   }
   
   curNode = &cpNode0;  
   while (curNode != NULL)
   {  if (confine2Plane[0] != 'Y') 
     for (i=0;i<curNode->n;i++) curNode->ary[i].c3 = randPM(initExtent,&idmy);  /* z coordinate */
      else
         for (i=0;i<curNode->n;i++) curNode->ary[i].c3 = 0.0;     
      curNode = curNode->nextnode;
   }
   /* now,  initialize velocities */
   curNode = &cpNode0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) curNode->ary[i].vx = randPM(initVelocity,&idmy);
      curNode = curNode->nextnode;
   }
   curNode = &cpNode0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) curNode->ary[i].vy = randPM(initVelocity,&idmy);
      curNode = curNode->nextnode;
   }
   curNode = &cpNode0; 
   while (curNode != NULL)
   {  if (confine2Plane[0] != 'Y')
     for (i=0;i<curNode->n;i++) curNode->ary[i].vz = randPM(initVelocity,&idmy);
      else 
         for (i=0;i<curNode->n;i++) curNode->ary[i].vz = 0.0;
      curNode = curNode->nextnode;
   }
}

void FixUpLargest000(  cmpNode *node0,double lFactor)
/* this gives the largest mass an unfair advantage */
{  cmpNode *curNode; int i;
  curNode = node0; while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) 
      if (curNode->ary[i].id == largestID)
      {  curNode->ary[i].m *= lFactor;
         curNode->ary[i].r *= sqrt(lFactor);  /* somewhat arbitrary */
         curNode->ary[i].c1 = 0.0;   /* init largest mass to be at the origin */
         curNode->ary[i].c2 = 0.0;
         curNode->ary[i].c3 = 0.0;
         
      }       
      curNode = curNode->nextnode;
   }
}



void displaceOrigin( double xDelta, double yDelta, double zDelta)
/* global: shifts centers of all bodies by xDelta,yDelta */
{ cmpNode *curNode; int i;
  curNode = &cpNode0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
     { curNode->ary[i].c1 += xDelta;
       curNode->ary[i].c2 += yDelta;
       curNode->ary[i].c3 += zDelta;
     }
      curNode = curNode->nextnode;
   }
} 


void OutputAllPoints(char *fStr)
/* global: outputs position,mass,radius,etc. of all points */
{
  cmpNode *curNode; int i;
  FILE *fH;
  fH = fopen(fStr,"w");
  curNode = &cpNode0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
     { fprintf(fH,"%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf\n",
       curNode->ary[i].c1,
       curNode->ary[i].c2,
       curNode->ary[i].c3,
       curNode->ary[i].id,
       curNode->ary[i].m,
       curNode->ary[i].r,
       curNode->ary[i].vx,
       curNode->ary[i].vy,
       curNode->ary[i].vz);
     }
      curNode = curNode->nextnode;
   }
   fclose(fH);
} 


/* used by below, possibly unneeded 
void clearGs(  cmpNode *node0)
{ cmpNode *curNode; int i;
  curNode = node0;  while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++) curNode->ary[i].gx = curNode->ary[i].gy = curNode->ary[i].gz 0.0;     
      curNode = curNode->nextnode;
   }
}
*/

void calcNewtonianGM(cmpNode *node0, double smallestRFac)
{ /* basic Newtonian gravity attraction, calc, where we retain for each body,
   the individual attractions (i.e. postpone summation for total attraction) */
   cmpNode *curNode; int i;
   double dx,dy,dz,dw,fmag,rsqrd,smallestRS;
   /* clearGs(node0); */
   curNode = node0; 
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)
       { /* dx = x0 - curNode->ary[i].c1;
         dy = y0 - curNode->ary[i].c2;
         dz = z0 - curNode->ary[i].c3; */
         dx = curNode->ary[i].c1; 
         dy = curNode->ary[i].c2;
         dz = curNode->ary[i].c3;
         rsqrd =  sqr(dx) + sqr(dy) + sqr(dz);
         
         if (rsqrd != 0.0)
         {  /* limit how close r exists */
           smallestRS = sqr(curNode->ary[i].r * smallestRFac);
           if (rsqrd < smallestRS) rsqrd = smallestRS;
           /* ubiquitous  1 / r^3/2 formula */
           fmag = curNode->ary[i].m / (rsqrd * sqrt(rsqrd));
               curNode->ary[i].gx = fmag * dx;
               curNode->ary[i].gy = fmag * dy;
               curNode->ary[i].gz = fmag * dz;
             }
             else
             { curNode->ary[i].gx = 0.0;   /* no self attraction ... rigid bodies */
               curNode->ary[i].gy = 0.0;
               curNode->ary[i].gz = 0.0;
             }
       }
                     
      curNode = curNode->nextnode;
   }
}

#define NEWTONS_G 6.67384e-11

  
#define MAX_DOUBLE 1.79e308 



void applyGravity( cmpNode *n0, int i0, double dt)
/* apply the effect of the gravitation of all bodies on body i0 (in buffer n0) ... 
   so we are going to bump each object at a time during the simulation, rather 
   than trying to effect infinitesimal change to all bodies at once */
{  cmpNode *curNode; int i;
  double gTx,gTy,gTz,ampF; 
  double xN1,xN2,xN3;
  double lBarrier2;
  gTx = gTy = gTz = 0.0;
   curNode = &cpNode0; /* start of buffer */
   while (curNode != NULL)
   {  for (i=0;i<curNode->n;i++)  /* if (curNode->ary[i].id != largestID) */
     {{ gTx += curNode->ary[i].gx ; 
        gTy += curNode->ary[i].gy ; 
        gTz += curNode->ary[i].gz ;
       } 
     }
      curNode = curNode->nextnode;
   }
   /* multiply by Newtons G */
   gTx *= NEWTONS_G;
   gTy *= NEWTONS_G;
   gTz *= NEWTONS_G;
   /*  printf("%d:Newton=%f %f ... \n ",i0,gTx,gTy); */
   xN1 = n0->ary[i0].c1;
   xN2 = n0->ary[i0].c2;
   xN3 = n0->ary[i0].c3;
   xN1 += dt * n0->ary[i0].vx;  /* v nought dt */
   xN2 += dt * n0->ary[i0].vy;
   xN3 += dt * n0->ary[i0].vz;
   xN1 += 0.5 * gTx * sqr(dt);
   xN2 += 0.5 * gTy * sqr(dt);
   xN3 += 0.5 * gTz * sqr(dt);
   /* bump the position of this body */
   n0->ary[i0].c1 = xN1;
   n0->ary[i0].c2 = xN2;
   n0->ary[i0].c3 = xN3;
   /* now, bump velocities */
   n0->ary[i0].vx += gTx * dt;
   n0->ary[i0].vy += gTy * dt;
   n0->ary[i0].vz += gTz * dt;
   
} 

void initCollG(  void )
{ cmpNode *curNode; int i; 
  curNode = &cpNode0; /* start of buffer */
  while (curNode != NULL)
  {  for (i=0;i<curNode->n;i++)  
       {  curNode->ary[i].cgx =0; /* clear collateral gravity attribute */
          curNode->ary[i].cgy =0; 
          curNode->ary[i].cgz =0;
       } 
     curNode = curNode->nextnode;
   }
} 


void Calc_CollF( cmpNode *epiNode, int epi )

/* epiNode and epi give the observation point (center body) */

{ cmpNode *curNode,*curNodeJ; int i,j;
  double dx,dy,dz,dw,fmag,rsqrd,smallestRS;
  double rAv1,rAv2,rAv3;
  double FC1,FC2,FC3;
  FC1 = 0.0; FC2 = 0.0; FC3 = 0.0;
  curNode = &cpNode0; while (curNode != NULL)
   { for (i=0;i<curNode->n;i++) if (curNode->ary[i].id != obsID)
     { curNodeJ = &cpNode0;
       while (curNodeJ != NULL)
       { /* for (j=0;j<curNodeJ->n;j++) */
    /* for now, since we only have one comp node working, just use j=i+1 */
        for (j=i+1;j<curNodeJ->n;j++)   /* temporary */  
         if ( (curNode->ary[i].id != curNodeJ->ary[j].id) && 
        (curNodeJ->ary[j].id != obsID))
         {  /* ok, end bodies are i and j, and the center body is at origin */
            rAv1 =  (curNode->ary[i].c1 +  curNodeJ->ary[j].c1) / 2.0;
              rAv2 =  (curNode->ary[i].c2 +  curNodeJ->ary[j].c2) / 2.0;
            rAv3 =  (curNode->ary[i].c3 +  curNodeJ->ary[j].c3) / 2.0;
        
         rsqrd =  ( sqr(rAv1) + sqr(rAv2) + sqr(rAv3) ) ;
         if (rsqrd != 0.0)
         {  /* limit how close r exists */
           smallestRS = sqr(epiNode->ary[epi].r);  /* collateral gravity does not happen if the center body intercepts path */
           if (rsqrd < smallestRS) rsqrd = smallestRS;
           /* coll gravity: outside body masses! ,   1 / r^3/2 formula */
           fmag = ( curNode->ary[i].m * curNodeJ->ary[j].m) / (rsqrd * sqrt(rsqrd));
           /* apply sigma later */
               FC1 += fmag * rAv1;
               FC2 += fmag * rAv2;
               FC3 += fmag * rAv3;            
         }
       }  /* end for loop ... j */
       curNodeJ = curNodeJ->nextnode;
     }
     
     } /* end for loop ... i */
     curNode = curNode->nextnode;
   }
   epiNode->ary[epi].cgx = FC1;
   epiNode->ary[epi].cgy = FC2;
   epiNode->ary[epi].cgz = FC3;
   
     
}


void applyCollGravity(cmpNode * curNode, int i, double collFract, double dt )
/* apply the effect of the collateral grav to all bodies at once */
{   double gTx,gTy,gTz,xN1,xN2,xN3;
  double sigma,m;
  
       {  sigma = (collFract * NEWTONS_G);
          m = curNode->ary[i].m;  /* mass of self */
          gTx = curNode->ary[i].cgx; /* collateral was assign directly to body attribute */
          gTy = curNode->ary[i].cgy ; 
          gTz = curNode->ary[i].cgz ;
     gTx *= (sigma / m); gTy *= (sigma / m); gTz *= (sigma / m);
     /* printf("%d:coll=%f %f ... ",i,gTx,gTy); */
     xN1 = curNode->ary[i].c1;
     xN2 = curNode->ary[i].c2;
     xN3 = curNode->ary[i].c3;
     xN1 += dt * curNode->ary[i].vx;  /* v nought dt */
     xN2 += dt * curNode->ary[i].vy;
     xN3 += dt * curNode->ary[i].vz;
     xN1 += 0.5 * gTx * sqr(dt);
     xN2 += 0.5 * gTy * sqr(dt);
     xN3 += 0.5 * gTz * sqr(dt);
     curNode->ary[i].c1 = xN1;/* bump the position of this body */
     curNode->ary[i].c2 = xN2;
     curNode->ary[i].c3 = xN3;
     curNode->ary[i].vx += gTx * dt; /* now, bump velocities */
     curNode->ary[i].vy += gTy * dt;
     curNode->ary[i].vz += gTz * dt;
       } 

} 



/* ************** Execution   *********** */

int main ( int argc, char *argv[])   

{
  int iEpoch,epi;
  cmpNode *epNode;
  int nbTotal,totalEpochs,i;
  double dt,initExtent,avgDensity,sigDensity,smRadius,incRadius,rFacBump;
  double degSwath;
  double globalOffset1,globalOffset2,globalOffset3;
  double lFactor,rBarrierFraction;
  char outFileStr[256];
  int initRandI;
  FILE *fPtrEP;
  /* get number of bodies, dt, fAmp, min sizes, etc. */
  printf("  Collateral gravity simulation\n");
  printf("    by Nathan E. Frick, summer 2011 \n");  
   
  printf("Enter the total number of bodies (e.g. 3000):\n");
  scanf("%d",&nbTotal);
  printf("Enter time increment in seconds (e.g. 10):\n");
  scanf("%lf",&dt);
  printf("Enter extent of simulation:\n");
  scanf("%lf",&initExtent);
  printf(" Confine to x-y plane (Y/N)?\n");
  scanf("%s",confine2Plane);
  printf("Enter extent of velocity seed (e.g. 1 m/s):\n");
  scanf("%lf",&initVelocity);
  printf("Enter average density in g/cc:\n");
  scanf("%lf",&avgDensity);
  printf("Enter density variance in g/cc:\n");
  scanf("%lf",&sigDensity);
  printf("Enter smallest radius to create (in meters):\n");
  scanf("%lf",&smRadius);
  printf("Enter radius increment (in meters):\n");
  scanf("%lf",&incRadius);
  printf("Enter radius bump for lognormal seed (e.g. 1.02) :\n");
  scanf("%lf",&rFacBump);
  

  printf("Enter factor to multiply times by G for collateral grav (e.g. 0.0043):\n");
  printf("This will be divided by a million before using ...\n");
  scanf("%lf",&collFraction);
  collFraction /= 1000000.0;
  printf("Enter number of epochs for simulation:\n");
  scanf("%d",&totalEpochs);
  
  printf("Enter direction swath (e.g. 4.0):\n");
  scanf("%lf",&degSwath);

  printf("Enter output filename:\n");
  scanf("%s",outFileStr);
  printf("Pre-multiply largest body by factor (e.g. 100):\n");
  scanf("%lf",&lFactor);
  printf("Enter an integer to seed random number generator for positions:\n");
  scanf("%d",&initRandI);
  /* convert units of density from g/cc to kg/m^3 */
  avgDensity *= 1000.0;
  sigDensity *= 1000.0;
  /* ok, now, begin */
  if (allocMem4Arrays(nbTotal) == 1) { printf("failed to allocate all mem\n"); return 1; }
  printf("initializing body radius, mass, position, etc ...\n");
  buildBodiesN((initExtent * 1.5),avgDensity,sigDensity,smRadius,incRadius,rFacBump,initRandI);
  FixUpLargest000(&cpNode0,lFactor);  /* initial mass increase */
  CheckForCollionsMerge();
  DumpZeroBodies(&cpNode0);
  OutputAllPoints(outFileStr);
  rBarrierFraction = 1.5;   /* (fraction of radius) distance where largest body -> no gravity */
  for (iEpoch=0;iEpoch < totalEpochs;iEpoch++)
  { printf("---------------------- EPOCH %d ----------------------\n",iEpoch);
    /*   if (iEpoch == 10) { FixUpLargest(&cpNode0,(1.0/lFactor)); } */  
     /* these can collapse on largest ?? */  
    initCollG();  /* must recalc each epoch */        
    epNode = &cpNode0; while (epNode != NULL)
     {  for (epi=0;epi<epNode->n;epi++) 
      {  /* epi is the current point of observation */
     obsID = epNode->ary[epi].id;
         globalOffset1 = epNode->ary[epi].c1;      
         globalOffset2 = epNode->ary[epi].c2;
         globalOffset3 = epNode->ary[epi].c3;
         displaceOrigin(-globalOffset1,-globalOffset2,-globalOffset3);
         /* calcBDirections(&cpNode0); */ 
         /* printf("global offset is %lf,%lf,%lf for observation point %d\n",
            globalOffset1,globalOffset2,globalOffset3,epNode->ary[epi].id); */
         /* printf("Calc of Newtonian gravity \n"); */
         calcNewtonianGM(&cpNode0,rBarrierFraction);  /* basic Newtonian gravity for this point */
         /* printf("Calc of collG\n"); */
         Calc_CollF(epNode,epi); 
         applyCollGravity(epNode,epi,collFraction,dt);
           /* printf("\n");
            for (i=0;i<cpNode0.n;i++) printf("for body %d: ampF=%lf \n",i,cpNode0.ary[i].ampF);  */       
         /* printf("Applying grav acceleration to point %d\r",epNode->ary[epi].id); */
         applyGravity(epNode,epi,dt);
         /* restore coordinates to original position */
         displaceOrigin(globalOffset1,globalOffset2,globalOffset3);
         /* PrintCollStats(); */
      }            
        epNode = epNode->nextnode;
        
     }
     OutputAllPoints(outFileStr);
     printf("epoch number %d is finished\n",iEpoch);
     fPtrEP = fopen("EPOCHNOW.TXT","w");
     fprintf(fPtrEP,"epoch number %d is finished\n",iEpoch);
     fclose(fPtrEP);
  /*  CheckForCollionsMerge();
     DumpZeroBodies(&cpNode0); */
     DumpFarBodies(&cpNode0,(initExtent * 1000.0));
  
  }


  return 0;
}







#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(int *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

float gasdev(int *idum)
/* returns a normally distrbuted deviate with zero mean
     and unit variance */
{
  static int iset=0;
  static float gset;
  float fac,r,v1,v2;
  float ran1();

  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0; /* pick two numbers from -1 to 1 */
      v2=2.0*ran1(idum)-1.0;
      r=v1*v1 + v2*v2;
    } while (r >= 1.0);  /* if they are not in the unity circle, try again */
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;  /* save for later */
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}



void initRand(int *idummy )

{ int idm;
  idm=-1;
  ran1(&idm);
  *idummy = idm;
}
  
float randPM (float maxOut,int *idummy)
/* plus maxOut to minus maxout */
{ return ((2*maxOut)*ran1(idummy)-maxOut);
}

float randP (float maxOut,int *idummy)
{ return (maxOut*ran1(idummy));
}

float gaussianDeviate (float mean,float variance, int *idummy)
{  return(mean + variance*gasdev(idummy));
}

