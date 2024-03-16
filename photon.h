/*c file photon_common.f*/
/*c*/
#ifndef PHOTONFLOATTYPE
#define PHOTONFLOATTYPE double
#endif
typedef PHOTONFLOATTYPE PHOTONFLOAT;
#define      npp0      5000
#define NPP npp0
#define ntasks0    1       
#define NTASKS ntasks0
#define ns0     37          
#define NSM ns0

#define   nseq0           1   

#define   tau0      0.100     
#define rhod0     0.300    
#define rhos0   0.500      

#define   modlus31  2147483647
#define epsdet0   1.e-10   
#define pi0    3.1415926536

#define   lensin0   90        
#define lencos0   4*lensin0
#define LENSIN lensin0
#define LENCOS lencos0
#define   nppm      npp0-1    
#define ntasksm   ntasks0-1 
#define nsm       ns0-1    

#define   nseqm     nseq0-1   
#define lensinm   lensin0-1 
#define lencosm   lencos0-1

/* flags for printp routine */
#define MAIN0 0
#define MAIN1 1
#define MAINF 0xf

/*c*/
/*c*/
/* not used on sun
      common/readonly/
  itask(0:1,0:ntasksm), jtask(0:ntasksm), time ,tarray(2)
 */
/*c*/
struct readonly
  {
    PHOTONFLOAT 
        x1      [NSM      ] ,x2      [NSM      ] ,
        y1      [NSM      ] ,y2      [NSM      ] ,
        delx    [NSM      ] ,dely    [NSM      ] ,
        enx     [NSM      ] ,eny     [NSM      ] ,
        sqln    [NSM      ] ,rhs     [NSM      ] ,
        spec1   [NSM      ] ,spec2   [NSM      ] ,
        sintab  [LENSIN  ] ,costab  [LENCOS  ] ;
  int
        itask[2][NTASKS], jtask[NTASKS];
  };
/*c*/
   /* ALLOCATE ONE FOR EACH TASK */
struct writeable
  {
    int iseedphi, iseedth, iseedrn, kevent, seconds;
    int itaske[NSM][NSM], itaskt[NSM][NSM];
    
  };

typedef struct writeable WRITEABLE;
typedef struct readonly READONLY;

