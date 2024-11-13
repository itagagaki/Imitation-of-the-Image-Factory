#define PLANE		1
#define	_2nd_ORDER	2

typedef struct _coord {
  double r11, r12, r13,
         r21, r22, r23,
         r31, r32, r33;
  double x0, y0, z0;
} COORD;

typedef struct _attr {
  int atno;
  bool tranflag;
  double diff_red, diff_grn, diff_blu;
  double spec_red, spec_grn, spec_blu;
  double tran_red, tran_grn, tran_blu;
  double lumi_red, lumi_grn, lumi_blu;
  double hilit, hilit_diff;
  double refract;
  int texmap;
} ATTR;

typedef struct _shape {
  int shno;
  int func;
  double a, b, c, d, e, f, g, h, i, j;
  vector position;
} SHAPE;

typedef struct and_term {
  bool positive;
  SHAPE *shape;
  ATTR *attribute;
  struct and_term *next_factor;
} ANDFACT;

typedef struct _orlist {
  ANDFACT *and_root;
  struct _orlist *next_term;
} ORLIST;

typedef struct _object {
  ORLIST *or_root;
  COORD coordinate;
} OBJECT;

typedef struct _cluster {
  ORLIST *or_root;
  double putx, puty, putz;
  struct _cluster *next_cluster;
} CLUSTER;

typedef struct _list {
  char *ptr;
  struct _list *next;
} LIST;
