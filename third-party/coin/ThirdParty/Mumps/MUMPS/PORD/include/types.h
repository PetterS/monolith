/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: types.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99sep14
/
/ This file contains the fundamental data structures
/
******************************************************************************/

typedef double FLOAT;
typedef int    options_t;
typedef FLOAT  timings_t;

/*****************************************************************************
Graph object
******************************************************************************/
typedef struct _graph {
  int   nvtx;
  int   nedges;
  int   type;
  int   totvwght;
  int   *xadj;
  int   *adjncy;
  int   *vwght;
} graph_t;

/*****************************************************************************
Graph bisection object
******************************************************************************/
typedef struct _gbisect {
  graph_t *G;
  int     *color;
  int     cwght[3];
} gbisect_t;

/*****************************************************************************
Domain decomposition object
******************************************************************************/
typedef struct _domdec {
  graph_t *G;
  int     ndom;
  int     domwght;
  int     *vtype;
  int     *color;
  int     cwght[3];
  int     *map;
  struct _domdec *prev, *next;
} domdec_t;

/*****************************************************************************
Bipartite graph object
******************************************************************************/
typedef struct _gbipart {
  graph_t *G;
  int     nX;
  int     nY;
} gbipart_t;

/*****************************************************************************
Recursive nested dissection object
******************************************************************************/
typedef struct _nestdiss {
  graph_t *G;
  int     *map;
  int     depth;
  int     nvint;
  int     *intvertex;
  int     *intcolor;
  int     cwght[3];
  struct _nestdiss *parent, *childB, *childW;
} nestdiss_t;

/*****************************************************************************
Multisector object
******************************************************************************/
typedef struct _multisector {
  graph_t *G;
  int     *stage;
  int     nstages;
  int     nnodes;
  int     totmswght;
} multisector_t;

/*****************************************************************************
Elimination graph object
******************************************************************************/
typedef struct _gelim {
  graph_t *G;
  int     maxedges;
  int     *len;
  int     *elen;
  int     *parent;
  int     *degree;
  int     *score;
} gelim_t;

/*****************************************************************************
Bucket structure object
******************************************************************************/
typedef struct _bucket {
  int   maxbin, maxitem;
  int   offset;
  int   nobj;
  int   minbin;
  int   *bin;
  int   *next;
  int   *last;
  int   *key;
} bucket_t;

/*****************************************************************************
Minimum priority object
******************************************************************************/
typedef struct _stageinfo stageinfo_t;
typedef struct _minprior {
  gelim_t       *Gelim;
  multisector_t *ms;
  bucket_t      *bucket;
  stageinfo_t   *stageinfo;
  int           *reachset;
  int           nreach;
  int           *auxaux;
  int           *auxbin;
  int           *auxtmp;
  int           flag;
} minprior_t;
struct _stageinfo {
  int   nstep;
  int   welim;
  int   nzf;
  FLOAT ops;
};

/*****************************************************************************
Elimination tree object
******************************************************************************/
typedef struct _elimtree {
  int   nvtx;
  int   nfronts;
  int   root;
  int   *ncolfactor;
  int   *ncolupdate;
  int   *parent;
  int   *firstchild;
  int   *silbings;
  int   *vtx2front;
} elimtree_t;

/*****************************************************************************
Input matrix object
******************************************************************************/
typedef struct _inputMtx {
  int   neqs;
  int   nelem;
  FLOAT *diag;
  FLOAT *nza;
  int   *xnza;
  int   *nzasub;
} inputMtx_t;

/*****************************************************************************
Dense matrix object
******************************************************************************/
typedef struct _workspace workspace_t;
typedef struct _denseMtx {
  workspace_t *ws;
  int         front;
  int         owned;
  int         ncol;
  int         nrow;
  int         nelem;
  int         nfloats;
  int         *colind;
  int         *rowind;
  int         *collen;
  FLOAT       *entries;
  FLOAT       *mem;
  struct _denseMtx *prevMtx, *nextMtx;
} denseMtx_t;
struct _workspace {
  FLOAT      *mem;
  int        size;
  int        maxsize;
  int        incr;
  denseMtx_t *lastMtx;
};

/*****************************************************************************
Compressed subscript structure object
******************************************************************************/
typedef struct _css {
  int   neqs;
  int   nind;
  int   owned;
  int   *xnzl;
  int   *nzlsub;
  int   *xnzlsub;
} css_t;

/*****************************************************************************
Front subscript object
******************************************************************************/
typedef struct _frontsub {
  elimtree_t *PTP;
  int        nind;
  int        *xnzf;
  int        *nzfsub;
} frontsub_t;

/*****************************************************************************
Factor matrix object
******************************************************************************/
typedef struct _factorMtx {
  int        nelem;
  int        *perm;
  FLOAT      *nzl;
  css_t      *css;
  frontsub_t *frontsub;
} factorMtx_t;

/*****************************************************************************
Mapping object
******************************************************************************/
typedef struct _groupinfo groupinfo_t;
typedef struct {
  elimtree_t  *T;
  int         dimQ;
  int         maxgroup;
  int         *front2group;
  groupinfo_t *groupinfo;
} mapping_t;
struct _groupinfo {
  FLOAT ops;
  int   nprocs;
  int   nfronts;
};

/*****************************************************************************
Topology object
******************************************************************************/
typedef struct {
  int      nprocs;
  int      mygridId;
  int      dimX;
  int      dimY;
  int      myQId;
  int      dimQ;
  int      *cube2grid;
#ifdef PARIX
  LinkCB_t **link;
#endif
#ifdef MPI
  MPI_Comm   comm;
  MPI_Status status;
#endif
} topology_t;

/*****************************************************************************
Communication buffer object
******************************************************************************/
typedef struct {
  char   *data;
  size_t len;
  size_t maxlen;
} buffer_t;

/*****************************************************************************
Bit mask object
******************************************************************************/
typedef struct {
  int   dimQ;
  int   maxgroup;
  int   mygroupId;
  int   offset;
  int   *group;
  int   *colbits, *colmask;
  int   *rowbits, *rowmask;
} mask_t;

