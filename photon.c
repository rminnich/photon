/*c*/
/*c*/
/*c                  *********************/
/*c                  *********************/
/*c                  **                ***/
/*c                  **  MAIN PROGRAM  ***/
/*c                  **                ***/
/*c                  *********************/
/*c                  *********************/
/*c*/
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "photon.h"

/* DO WE NEED IFIX FOR REAL? */
#define ifix(f) (int) (f)

#ifdef MAIN
int
itotale [NSM][NSM] ,itotalt [NSM][NSM] ,
irowsume[NSM      ] ,irowsumt[NSM      ];
main(argc, argv)
int argc;
char *argv[];
{
  static struct readonly readonly[1];
  static struct writeable writeable[16];
  struct writeable *lw;
  int itsk;
  int nevents;
  PHOTONFLOAT time;
  int li, le;
  printf("RUNNING SUN VERSION.\n");

  /* now separate to remove assumptions about contiguous data */
  initreadonly(readonly);
  for(itsk = 0; itsk < ntasks0; itsk++)
    initwriteable(itsk, &writeable[itsk]);

  prntp(MAIN0, 0, readonly, writeable);
  for(itsk = 0; itsk < NTASKS; itsk++)
    taskcode(readonly->jtask[itsk], readonly, &writeable[itsk]);
  time = nevents = 0;
  for(itsk = 0, lw=writeable; itsk < NTASKS; itsk++, lw++)
    {
              nevents += lw->kevent;
              time += lw->seconds;
              for(le = 0; le < NSM; le++)
                for(li = 0; li < NSM; li++)
                  {
                    itotale[le][li] =
                      itotale[le][li] + lw->itaske[le][li];
                    itotalt[le][li] =
                      itotalt[le][li] + lw->itaskt[le][li];
                  }
    }

  done(nevents, time, readonly);
}
#else
extern int
itotale [NSM][NSM] ,itotalt [NSM][NSM] ,
irowsume[NSM      ] ,irowsumt[NSM      ];
#endif
/*c            *********************************/
/*c            *********************************/
/*c            **                            ***/
/*c            **  CODE FOR A SINGLE THREAD  ***/
/*c            **                            ***/
/*c            *********************************/
/*c            *********************************/
/*c*/
/*c*/
taskcode(jtsk, readonly, writeable)
int jtsk;
struct readonly *readonly;
struct writeable *writeable;
{
  PHOTONFLOAT dxcrnt = 0.0, dycrnt = 0.0, toffset = 0.0, xe0 = 0.0, ye0 = 0.0, 
  exloc = 0.0, eyloc = 0.0, ex = 0.0, ey = 0.0, delxl = 0.0, delyl = 0.0, 
  rhsl = 0.0, absdt = 0.0, scldenom = 0.0, scalephi = 0.0, scaletha = 0.0, 
  fitsk = 0.0, fntasks = 0.0, fnplfinv = 0.0, xe = 0.0, ye = 0.0,
  rhse = 0.0, det = 0.0, dtinv = 0.0, xi = 0.0, yi = 0.0, x1l = 0.0, y1l = 0.0, 
  y1 = 0.0, x2l = 0.0, x2 = 0.0, y2l = 0.0, y2 = 0.0, sqlnl = 0.0, sqln = 0.0, 
  ssq = 0.0, rani = 0.0;
  PHOTONFLOAT exold = 0.0, eyold = 0.0, spec1c2 = 0.0, spec2c2 = 0.0, eny = 0.0, 
  enx = 0.0;
  int itsk= 0, le0= 0, ixphi= 0, ixth= 0, iphi= 0, ith= 0, idead= 0, li= 0, 
   l= 0, i= 0, le= 0, ix = 0;
  int timestart = 0, timeend = 0;
#define sun
#ifdef sun
  struct rusage rusage[1];

  if (getrusage(RUSAGE_SELF, rusage) < 0)
    perror("rusage fails");
  timestart = rusage[0].ru_utime.tv_sec;
#endif
  itsk = jtsk;
  printf("in taskcode - itsk= %d\n",itsk);
  printf("readonly is 0x%x, writeable is 0x%x\n", readonly, writeable);

  scldenom = 1.0/((PHOTONFLOAT)(modlus31)-2.0);
  scalephi  = (PHOTONFLOAT)(lencos0)*scldenom;
  scaletha  = (PHOTONFLOAT)(lensin0)*scldenom;
  /*c*/
  fitsk    = (PHOTONFLOAT)(itsk);
  fntasks  = (PHOTONFLOAT)(ntasks0);
  fnplfinv = 1./(PHOTONFLOAT)(npp0);
  /*c*/
  for(le0 = 0; le0 < NSM; le0++)
    {
      dxcrnt   = readonly->delx[le0]*fnplfinv;
      dycrnt   = readonly->dely[le0]*fnplfinv;
      toffset  = (fitsk + .5)/fntasks;
      xe0      = readonly->x1[le0] + (toffset - 1.0)*dxcrnt;
      ye0      = readonly->y1[le0] + (toffset - 1.0)*dycrnt;
      /*c*/
      for(i = 0; i < NPP; i++)
	{
	  xe0      = xe0 + dxcrnt;
	  ye0      = ye0 + dycrnt;
	  xe       = xe0;
	  ye       = ye0;
	  le       = le0;
	  /*c*/
	  /*c                        EMIT*/
	  /*c*/
	  writeable->iseedphi= iran31(writeable->iseedphi);
	  writeable->iseedth = iran31(writeable->iseedth );
	  ixphi         = writeable->iseedphi - 1;
	  ixth          = writeable->iseedth  - 1;
	  iphi          = ifix((PHOTONFLOAT)(ixphi)*scalephi);
	  ith           = ifix((PHOTONFLOAT)(ixth )*scaletha);
	  exloc         = readonly->sintab[ith]*readonly->costab[iphi];
	  eyloc         = readonly->costab[ith];
	  ex            = exloc*readonly->eny[le] + eyloc*readonly->enx[le];
	  ey            = eyloc*readonly->eny[le] - exloc*readonly->enx[le];
	  idead         = 0;
	  /*c*/
	  /*c*/
	  /*c*/
	  /*c                        INTSEC*/
	  /*c*/
	  /*c*/
	  /*c         *    main loop over surfaces    **/
	  /*c   * for now this only works for convex surfaces **/
	  /*c*/
	  while (! idead)
	    {
	      writeable->kevent  = writeable->kevent + 1;
	      rhse          = ex*ye - ey*xe;
	      li = le;
	      for(l = 0; l < NSM; l++)
		{
		  if(l != le) 
		    {
		      delxl   = readonly->delx   [l];
		      delyl   = readonly->dely   [l];
		      rhsl    = readonly->rhs    [l];
		      /*c*/
		      /*c  compute intersection points*/
		      /*c*/
		      det  = ex*delyl - ey*delxl;
		      absdt= fabs(det);
		      if(absdt <= epsdet0) 
			det=  epsdet0;
		      dtinv= 1.0/det;
		      xi     = dtinv * (delxl*rhse - ex*rhsl);
		      yi     = dtinv * (delyl*rhse - ey*rhsl);
		      /*c*/
		      /*c  test for intersection between surface endpoints*/
		      /*c*/
		      x1l     = readonly->x1     [l];
		      y1l     = readonly->y1     [l];
		      x2l     = readonly->x2     [l];
		      y2l     = readonly->y2     [l];
		      sqlnl   = readonly->sqln   [l];
		      ssq  = (xi - x1l)*(xi - x1l) + (xi - x2l)*(xi - x2l)
			+ (yi - y1l)*(yi - y1l) + (yi - y2l)*(yi - y2l);
		      if(ssq <= sqlnl) 
			{
			  li= l;
			  break;
			}
		    }
		}

/*c*/
/*c         *    end of loop over surfaces    **/
/*c*/
/*c*/
/*c*/
/*c*/
/*c                            ABSREF*/
/*c*/
/*c ****************************************************************/
/*c *                                                             **/
/*c *                 absref disposition:                         **/
/*c *                                                             **/
/*c *            1.0 -                                            **/
/*c *                | ===> absorption in this interval (e++)     **/
/*c *           rhos -                                            **/
/*c *                | ===> specular reflection in this interval  **/
/*c *  rani     rhod -                                            **/
/*c *                | ===> diffuse reflection in this interval   **/
/*c *            tau -                                            **/
/*c *                | ===> transmission in this interval (t++)   **/
/*c *            0.0 -                                            **/
/*c *                                                             **/
/*c ****************************************************************/
/*c*/
	      if(li != le)
		{
		  /*c*/
		  /*c  intersection found*/
		  /*c*/
		  writeable->iseedrn = iran31(writeable->iseedrn );
		  ix            = writeable->iseedrn  - 1;
		  rani          = (PHOTONFLOAT)(ix)*scldenom;
		  /*c*/
		  /*c  test for absorbed photon*/
		  /*c*/
		  if(rani >= rhos0)
		    {
		      writeable->itaske[le0][li] = 
			writeable->itaske[le0][li] + 1;
		      idead = 1;
		    }
		  else
		    /*c*/
		    /*c  test for specularly reflected photon, and reemit*/
		    /*c*/
		    if(rani >= rhod0)
		      {
			exold  = ex;
			eyold  = ey;
			xe     = xi;
			ye     = yi;
			spec1c2= readonly->spec1[li];
			spec2c2= readonly->spec2[li];
			ex     = spec1c2*exold + spec2c2*eyold;
			ey     = spec2c2*exold - spec1c2*eyold;
			le     = li;
		      }
		  else
		    /*c*/
		    /*c  test for diffusely reflected photon*/
		    /*c*/
		    if (rani >= tau0)
		      {
			xe            = xi;
			ye            = yi;
			le            = li;
			writeable->iseedphi = iran31(writeable->iseedphi);
			writeable->iseedth = iran31(writeable->iseedth);
			ixphi         = writeable->iseedphi - 1;
			ixth          = writeable->iseedth  - 1;
			iphi          = ifix((PHOTONFLOAT)(ixphi)*scalephi);
			ith           = ifix((PHOTONFLOAT)(ixth )*scaletha);
			exloc         = readonly->sintab[ith]*
			  readonly->costab[iphi];
			eyloc         = readonly->costab[ith];
			ex            = exloc*readonly->eny[le] + 
			  eyloc*readonly->enx[le];
			ey            = eyloc*readonly->eny[le] - 
			  exloc*readonly->enx[le];
		      }
		  else
		    {
		      /*c*/
		      /*c  transmitted photon identified*/
		      /*c*/
		      writeable->itaskt[le0][li]= 
			writeable->itaskt[le0][li] + 1;
		      idead = 1;
		    }
		}
	      else
		{
		  /*c*/
		  /*c  no intersection found --- photon is lost*/
		  /*c*/
		  idead = 1;
		}
/*c*/
	    }
	}
    }
#ifdef sun
  getrusage(RUSAGE_SELF, rusage);
  timeend = rusage[0].ru_utime.tv_sec;
#endif
  printf("total task time start 0x%x end 0x%d diff 0x%d\n", 
          timestart, timeend, timeend-timestart);
  writeable->seconds = timeend-timestart;
  fflush(stdout);
/**/
  prntp(MAIN1, itsk, readonly, writeable);
/* */

}
/*c*/
/*c*/
/*c*/
/*c*/
/*c*/
initreadonly(readonly)
struct readonly *readonly;
  {
    PHOTONFLOAT dxl= 0.0, dyl= 0.0, scale = 0.0;
    int itsk= 0, j= 0, i= 0, l = 0;
    extern void getngon();
    printf("initreadonly with 0x%x\n", readonly);
    fflush(stdout);
    for(itsk = 0; itsk < NTASKS; itsk++)   
      {
         readonly->itask[0][itsk] = 2;
         readonly->jtask  [itsk] = itsk;
      }

    getngon(readonly->x1,readonly->y1,readonly->x2,readonly->y2);
/*c*/
    for(j = 0; j < NSM; j++)
      for(i = 0; i < NSM; i++)
	{
	  itotale[i][j     ]  = 0;
	  itotalt[i][j     ]  = 0;
	}
    
    /*c*/
    for(l = 0; l < NSM; l++)
      {
         dxl        = readonly->x2[l] - readonly->x1[l];
         dyl        = readonly->y2[l] - readonly->y1[l];
         readonly->delx    [l]= dxl;
         readonly->dely    [l]= dyl;
         readonly->sqln    [l]= dxl*dxl + dyl*dyl;
         scale      = 1.0/sqrt(readonly->sqln[l]);
         readonly->enx     [l]=  dyl*scale;
         readonly->eny     [l]= -dxl*scale;
         readonly->rhs     [l]= dxl*readonly->y1[l] - dyl*readonly->x1[l];
         readonly->spec1   [l]= readonly->eny[l]*readonly->eny[l] - 
	   readonly->enx[l]*readonly->enx[l];
         readonly->spec2   [l]= -2.0*readonly->eny[l]*readonly->enx[l];
         irowsume[l]= 0;
         irowsumt[l]= 0;
       }
/*c*/
    for(i = 0; i < LENSIN; i++)
      readonly->sintab[i]= sin(((PHOTONFLOAT)(i) + .5)*pi0/(2.0*(PHOTONFLOAT)(lensin0)));
/*c*/
    for(i = 0; i < LENCOS; i++)
      readonly->costab[i]= cos(((PHOTONFLOAT)(i) + .5)*pi0*2.0/(PHOTONFLOAT)(lencos0));

  }
initwriteable(itsk, wr)
int itsk;
struct writeable *wr;
  {
    register struct writeable *lw = wr;
    int j, i;
    printf("initwriteable task %d addr 0x%x\n", itsk, wr);

    lw->iseedphi = 3*(itsk) + 1;
    lw->iseedth  = 3*(itsk) + 2;
    lw->iseedrn = 3*(itsk) + 3;
    
    lw->kevent = 0;
    for(j = 0; j < NSM; j++)
      for(i = 0; i < NSM; i++)
	{
	  lw->itaske [i][j]  = 0;
	  lw->itaskt [i][j]  = 0;
	}

  }
/*c*/
/*c*/
/*c*/
/*c*/
/*c*/
void
getngon(x1,y1,x2,y2)
PHOTONFLOAT 
x1[],y1[],x2[],y2[];
{
  int ns = 0, l = 0;
  PHOTONFLOAT scale = 0., dr = 0., thl = 0.;
  /*c*/
  /*c  construct outer polygon (counter clockwise)*/
  /*c*/
  ns   = nsm+1;
  scale= 1.;
  dr   = -2.*3.14159265358/(PHOTONFLOAT)(ns);
  /*c*/
  x1[0]= scale*1.;
  y1[0]= scale*0.;
  for(l = 1; l < NSM; l++)
    {
      thl= (PHOTONFLOAT)(l)*dr;
      x1[l]= scale*cos(thl);
      y1[l]= scale*sin(thl);
      x2[l-1]= x1[l];
      y2[l-1]= y1[l];
    }
  x2[nsm]= x1[0];
  y2[nsm]= y1[0];
}
/*c*/
/*c*/
/*c*/
/*c*/
/*c*/
prntp(iflag, istop, readonly, wr)
int iflag, istop;
struct readonly *readonly;
struct writeable *wr;
{
  int nevents = 0, iesum = 0, itotal = 0, itsk = 0, le = 0, li = 0;

  struct writeable *lw;
  PHOTONFLOAT evperpht= 0, perphtn= 0, time= 0, itsum = 0;
  if (iflag == MAIN0)
    {
      
         printf("        i n i t i a l   d a t a \n");
         printf(" nseq0         = %d\n",nseq0);
         printf(" ntasks0       = %d\n"   ,ntasks0);
         printf(" ns0           = %d\n"   ,ns0);
         printf(" npp0          = %d\n"   ,npp0);
         printf(" epsdet0       = %f\n",epsdet0);
         printf(" tau0          = %f\n", tau0);
         printf(" rhod0         = %f\n",rhod0);
         printf(" rhos0         = %f\n",rhos0);
         printf(" lensin0       = %d\n",lensin0);
         printf(" lencos0       = %d\n",lencos0);
         printf("e n d   o f   i n i t i a l   d a t a\n");
	fflush(stdout);
         if(istop)
           return;
       }
  

  /*c*/
  /*c*/
  if (iflag == MAIN1)
    {
      /*c*/
      /*c  using istop for task number ---*/
      /*c*/
      itsk = istop;
      lw = wr;
      nevents = lw->kevent;
      iesum   = 0;
      itsum   = 0;
      for(le = 0; le < NSM; le++)
        for(li = 0; li < NSM; li++)
          {
            iesum      = iesum       + lw->itaske[le][li];
            itsum      = itsum       + lw->itaskt[le][li];
          }
      itotal     = iesum + itsum;
      evperpht = (PHOTONFLOAT)(nevents)/(PHOTONFLOAT)(itotal);
      /*c*/
      printf(" *** final output -- itsk = %d\n",itsk);
      /*c*/
      if (nsm < 10)
        {
          printf(" itaske = ");
          for(le = 0; le < NSM; le++)
            for(li = 0; li < NSM; li++)
              printf("%d ",(lw->itaske[le][li]));
          printf("\n");
        }

      if ((nsm < 10) && (itsum > 0.5))
        {
          printf(" itaskt = ");
          for(le = 0; le < NSM; le++)
            for(li = 0; li < NSM; li++)
              printf("%d ",(lw->itaskt[le][li]));
          printf("\n");
        }


/*c*/
    printf("\n");
    printf(" number of events                = %d\n"   ,nevents);
    printf(" number of events/photon         = %f\n" ,evperpht);
    printf(" number of photons absorbed      = %d\n"   ,iesum);
    printf(" number of photons transmitted   = %d\n"   ,itsum);
    printf(" total number of photons counted = %d\n"   ,itotal);
    printf(" nevents = %d\n"   ,nevents);
    printf(" nev/phn = %f\n",(PHOTONFLOAT)(nevents)/(PHOTONFLOAT)(itotal));
    printf(" iesum   = %d\n"   ,iesum);
    printf(" itsum   = %d\n"   ,itsum);
    printf(" itotal  = %d\n"   ,itotal);
  }
/*c*/
/*c*/
}
done(nevents, time, readonly)
int nevents; 
PHOTONFLOAT time;
struct readonly *readonly;
{
  int iesum = 0, itotal = 0, itsk = 0, le = 0, li = 0;

  PHOTONFLOAT evperpht= 0, perphtn= 0, itsum = 0;


/*c*/
	  iesum   = 0;
	  itsum   = 0;
	  for(le = 0; le < NSM; le++)
	    {
	      irowsume[le]= 0;
	      irowsumt[le]= 0;
	      for(li = 0; li < NSM; li++)
		{
		  irowsume[le]= 
		    irowsume[le] + itotale[le][li];
		  irowsumt[le]= 
		    irowsumt[le] + itotalt[le][li];
		  iesum       = iesum        + itotale[le][li];
		  itsum       = itsum        + itotalt[le][li];
		}
	    }
	  itotal    = iesum + itsum;
	  evperpht  = (PHOTONFLOAT)(nevents)/(PHOTONFLOAT)(itotal);
	  perphtn   = time/(PHOTONFLOAT)(itotal);
	  /*c*/
	  printf(" *** final output ***\n");
	  /*c*/
	  if(nsm < 100)
	    {
	      printf(" irowsume = ");
	      for(le = 0; le < NSM; le++)
		printf(" %d ", irowsume[le]);
	      printf("\n");
	      if(itsum >= 0.5)
		{
		  printf(" irowsumt = ");
		  for(le = 0; le < NSM; le++)
		    printf(" %d ", irowsumt[le]);
		  printf("\n");
		  
		}
	    }
	  /*c*/
	  if (nsm < 10)
	    {
	      printf(" itotale = ");
	      for(le = 0; le < NSM; le++)
		for(li = 0; li < NSM; li++)
		  printf("%d ",(itotale[le][li]));
	      
	    }
	  if ((nsm < 10) && (itsum > 0.5))
	    {
	      printf(" itotalt = ");
	      for(le = 0; le < NSM; le++)
		for(li = 0; li < NSM; li++)
		  printf("%d ",(itotalt[le][li]));
	      
	    }
	  /*c*/
/*c*/
      printf("\n");
      printf(" number of events                = %d\n", nevents);
      printf(" number of events/photon         = %f\n", evperpht);
      printf(" number of photons absorbed      = %d\n", iesum);
      printf(" number of photons transmitted   = %f\n", itsum);
      printf(" total number of photons counted = %d\n", itotal);
      printf(" CPU time (sec)                  = %f\n", time);
      printf(" CPU time/photon (sec)           = %f\n", perphtn);
}

/*c*/
/*c*/
stop(c)
char *c;
{
  printf("%s\n", c);
	fflush(stdout);
  exit(1);
}
/*c*/
/*c*/
/*c*/
/*c*/
/*c*/
int
  iran31(oldseed)
int oldseed;
{
  /*c*/
  /*c returns a random integer in {1,2,...,2147483646 = 2^31 - 2}*/
  /*c -- see Park & Miller, Communications, ACM Nov, 1988*/
  /*c*/
  int hi,lo,seed;
#define a_right  16807
#define q_31     127773
#define   r_31     2836 
#define    m_31     2147483647
  /*c*/
  hi     = oldseed/q_31;
  lo     = oldseed - (hi*q_31);
  seed   = a_right*lo - r_31*hi;
  if(seed <- 0) 
    seed = seed + m_31;
  return seed;
  
}
