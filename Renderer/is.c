#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include "config.h"
#include "form.h"
#include "vector.h"
#include "object.h"
#include "light.h"
#include "intcode.h"
#include "input.h"


int reso_h = 256;
int reso_v = 240;
int cam_in_obj_flag = 0;
int obj_complex_flag = 0;

int
  hilit_sw = 1,
  shadow_sw = 1;

int trace_count = 20;
double min_energy = 0.02;
int super_level = 0;
double ambience = 0.0;

double min_expseed;

double pixel_size_x, pixel_size_y;
double pixel_size_xs, pixel_size_ys;
double screen_left_edge, screen_upper_edge;

char inpname[64], outname[64], maskout[64];
FILE *wfp;

int anti_depth;
int anti_thres;
int num_light;

char light_kind[MAXNL];
double light_x[MAXNL],light_y[MAXNL],light_z[MAXNL],
       light_vx[MAXNL],light_vy[MAXNL],light_vz[MAXNL],
       light_r[MAXNL], light_s[MAXNL],
       light_red[MAXNL],light_grn[MAXNL],light_blu[MAXNL];

double camera_x, camera_y, camera_z, focal_dist,
       ca11, ca12, ca13, ca21, ca22, ca23,
       ca31z, ca32z, ca33z;

CLUSTER *object_root, *last_object;
LIST *shl_root, *last_shl, *atl_root, *last_atl;


void warning(s)
char *s;
{
  fprintf(stderr, "warning: %s\n", s);
}

void init_tree()
{
  object_root = NULL;
  shl_root = atl_root = NULL;
  num_light = 0;
}

static void free_list(root)
LIST *root;
{
  LIST *list;
  char *freeptr;

  for (list = root; list; ) {
    free((char *)list->ptr);
    freeptr = (char *)list;
    list = list->next;
    free(freeptr);
  }
}

void clrpara()
{
  CLUSTER *objptr;
  ORLIST *orptr;
  ANDFACT *andptr;

  for (objptr = object_root; objptr; ) {
    char *freeptr;

    freeptr = (char *)objptr;
    for (orptr = objptr->or_root; orptr; ) {
      char *freeptr;

      freeptr = (char *)orptr;
      for (andptr = orptr->and_root; andptr; ) {
	char *freeptr;

	freeptr = (char *)andptr;
	andptr = andptr->next_factor;
	free(freeptr);
      }
      orptr = orptr->next_term;
      free(freeptr);
    }
    objptr = objptr->next_cluster;
    free(freeptr);
  }

  free_list(shl_root);
  free_list(atl_root);
  init_tree();
}

static char *trace_list(root, length)
LIST *root;
int length;
{
  register int count;
  register LIST *list;

  list = root;
  for (count = 1; count < length && list; count++, list = list->next) ;

  if (count != length) {
    char s[16];
    sprintf(s, "bad list %2d %2d", count, length);
    warning(s);
  }

  return list->ptr;
}

void setpara(fp)
FILE *fp;
{
  int ch;

  while ((ch = getc(fp)) != EOF)
    switch (ch)
      {
      case INTOBJECT:
	{
	  register int leafno;
	  register CLUSTER *object;

	  object = (CLUSTER *)malloc(sizeof(CLUSTER));
	  if (object_root) last_object->next_cluster = object;
	  else object_root = object;

	  object->next_cluster = NULL;
	  object->putx = getd(fp);
	  object->puty = getd(fp);
	  object->putz = getd(fp);
	  object->or_root = NULL;

	  while (getc(fp) == 0) {
	    ORLIST *term, *last_term;

	    term = (ORLIST *)malloc(sizeof(ORLIST));
	    if (object->or_root) last_term->next_term = term;
	    else object->or_root = term;

	    term->next_term = NULL;
	    term->and_root = NULL;

	    while (getc(fp) == 0) {
	      ANDFACT *factor, *last_factor;
	      LIST *list;

	      factor = (ANDFACT *)malloc(sizeof(ANDFACT));
	      if (term->and_root) last_factor->next_factor = factor;
	      else term->and_root = factor;
	      factor->next_factor = NULL;

	      if ((leafno = geti(fp)) == 0) {
		register SHAPE *sh;

		sh = (SHAPE *)malloc(sizeof(SHAPE));
		list = (LIST *)malloc(sizeof(LIST));
		if (shl_root) last_shl->next = list;
		else shl_root = list;

		list->ptr = (char *)sh;
		list->next = NULL;
		last_shl = list;

		sh->func = getc(fp);
		sh->a = getd(fp);
		sh->b = getd(fp);
		sh->c = getd(fp);
		sh->d = getd(fp);
		sh->e = getd(fp);
		sh->f = getd(fp);
		sh->g = getd(fp);
		sh->h = getd(fp);
		sh->i = getd(fp);
		sh->j = getd(fp);
		sh->position.x = getd(fp);
		sh->position.y = getd(fp);
		sh->position.z = getd(fp);
		factor->shape = sh;
	      }
	      else factor->shape = (SHAPE *)trace_list(shl_root, leafno);

	      if ((leafno = geti(fp)) == 0) {
		register ATTR *at;

		at = (ATTR *)malloc(sizeof(ATTR));
		list = (LIST *)malloc(sizeof(LIST));
		if (atl_root) last_atl->next = list;
		else atl_root = list;

		list->ptr = (char *)at;
		list->next = NULL;
		last_atl = list;

		at->tranflag = getc(fp);
		at->diff_red = getd(fp);
		at->diff_grn = getd(fp);
		at->diff_blu = getd(fp);
		at->spec_red = getd(fp);
		at->spec_grn = getd(fp);
		at->spec_blu = getd(fp);
		at->tran_red = getd(fp);
		at->tran_grn = getd(fp);
		at->tran_blu = getd(fp);
		at->lumi_red = getd(fp);
		at->lumi_grn = getd(fp);
		at->lumi_blu = getd(fp);
		at->hilit = getd(fp);
		at->hilit_diff = getd(fp);
		at->refract = getd(fp);
		at->texmap = getc(fp);
		factor->attribute = at;
	      }
	      else factor->attribute = (ATTR *)trace_list(atl_root, leafno);

	      factor->positive = getc(fp);
	      last_factor = factor;
	    }

	    last_term = term;
	  }

	  last_object = object;
	}
	break;

      case INTCAMERA:
	camera_x = getd(fp);
	camera_y = getd(fp);
	camera_z = getd(fp);
	focal_dist = getd(fp);
	ca11 = getd(fp);
	ca12 = getd(fp);
	ca13 = getd(fp);
	ca21 = getd(fp);
	ca22 = getd(fp);
	ca23 = getd(fp);
	ca31z = getd(fp);
	ca32z = getd(fp);
	ca33z = getd(fp);
	break;

      case INTLIGHT:
	if (num_light >= MAXNL) warning("too many lights");
	switch (light_kind[num_light] = getc(fp))
	  {
	  case PARA_LIT:
	    light_vx[num_light] = getd(fp);
	    light_vy[num_light] = getd(fp);
	    light_vz[num_light] = getd(fp);
	    break;

	  case POI_LIT:
	  case FPOI_LIT:
	    light_x[num_light] = getd(fp);
	    light_y[num_light] = getd(fp);
	    light_z[num_light] = getd(fp);
	    break;

	  case SPOT_LIT:
	    light_vx[num_light] = getd(fp);
	    light_vy[num_light] = getd(fp);
	    light_vz[num_light] = getd(fp);
	    light_x[num_light] = getd(fp);
	    light_y[num_light] = getd(fp);
	    light_z[num_light] = getd(fp);
	    light_r[num_light] = getd(fp);
	    light_s[num_light] = getd(fp);
	    break;

	  case AMBI_LIT:
	  default:
	    break;
	  }

	light_red[num_light] = getd(fp);
	light_grn[num_light] = getd(fp);
	light_blu[num_light++] = getd(fp);
	break;

      default:
	printf("%x\n", ch);
	warning("Undefined data");
      }
}

int atouie(s, termc, nextptr)
char *s, **nextptr;
int termc;
{
  register int val = 0;

  if (*s == termc) val = -1;
  else {
    while (isdigit(*s))
      if ((val = 10 * val + *s++ - '0') < 0) {
	val = -1;
	break;
      }

    if (*s != termc) val = -1;
  }
  *nextptr = ++s;

  return val;
}

void main(argc, argv)
int argc;
char *argv[];
{
  extern void cam_refract();
  extern void compute();

  int aspect_h, aspect_v;
  double shift_h, shift_v;
  double screen_size_x, screen_size_y;
  char iname[64];
  FILE *fp;



  min_expseed = log(TINY_VAL);

  init_tree();


/****************************
 *                          *
 *   COMMAND LINE ANALYZE   *
 *                          *
 ****************************/

  {
    bool outdes = FALSE;
    bool resdes = FALSE;
    bool aspdes = FALSE;
    bool trcdes = FALSE;
    bool mindes = FALSE;

    if (argc <= 1) {
      fprintf(stderr, "\nUsage:- ira [-flags] <filename>\n");
      fprintf(stderr, "flags:-\n");
      fprintf(stderr, "        a - specificate anti-aliasing level             -a <n>\n");
      fprintf(stderr, "        h - un-highlighting\n");
      fprintf(stderr, "        i - camera is in the object\n");
      fprintf(stderr, "        l - specificate minimum limit of energy of ray  -l <m>\n");
      fprintf(stderr, "        o - designate output                            -o <filename>\n");
      fprintf(stderr, "        r - specificate resolution                      -r <H>*<V>\n");
      fprintf(stderr, "        s - un-shadowing\n");
      fprintf(stderr, "        t - specificate maximum limit of trace level    -t <n>\n");
      fprintf(stderr, "        v - specificate aspect ratio                    -v <H>:<V>\n");
      exit(1);
    }

    while (--argc > 0 && **++argv == '-') {
      char *flagptr = *argv;

      while (*++flagptr) {
	switch (*flagptr)
	  {
	  case 'a':
	    if (--argc == 0) {
	      fprintf(stderr, "missing specification of anti-aliasing level\n");
	      exit(1);
	    }

	    {
	      char *descptr;

	      descptr = *++argv;
	      if ((super_level = atouie(descptr, '\0', &descptr)) < 0) {
		fprintf(stderr, "invalid anti-aliasing level description %s\n", *argv);
		exit(1);
	      }
	    }
	    break;

	  case 'h':
	    hilit_sw = 0;
	    break;

	  case 'i':
	    cam_in_obj_flag = 1;
	    break;

	  case 'c':
	    obj_complex_flag = 1;
	    break;

	  case 'l':
	    if (mindes) {
	      fprintf(stderr, "multiple specification of minimum limit of energy\n");
	      exit(1);
	    }
	    if (--argc == 0) {
	      fprintf(stderr, "missing specification of minimum limut of energy\n");
	      exit(1);
	    }
	    if ((min_energy = atof(*++argv)) < 0.0) {
	      fprintf(stderr, "invalid min_energy %s\n", *argv);
	      exit(1);
	    }
	    if (min_energy > 1.0) {
	      fprintf(stderr, "invalid min_energy %s\n", *argv);
	      exit(1);
	    }
	    mindes = TRUE;
	    break;

	  case 'o':
	    if(outdes) {
	      fprintf(stderr, "multiple output designation flag\n");
	      exit(1);
	    }
	    if (--argc == 0) {
	      fprintf(stderr, "missing description of output\n");
	      exit(1);
	    }
	    strcpy(outname, *++argv);
	    outdes = TRUE;
	    break;

	  case 'r':
	    if (resdes) {
	      fprintf(stderr, "multiple resolution spesification flag\n");
	      exit(1);
	    }
	    if (--argc == 0) {
	      fprintf(stderr, "missing specification of resolution\n");
	      exit(1);
	    }

	    {
	      char *descptr = *++argv;

	      if ((reso_h = atouie(descptr, '*', &descptr)) < 1) {
		fprintf(stderr, "invalid resolution description %s\n", *argv);
		exit(1);
	      }
	      if ((reso_v = atouie(descptr, '\0', &descptr)) < 1) {
		fprintf(stderr, "invalid resolution description %s\n", *argv);
		exit(1);
	      }
	    }
	    resdes = TRUE;
	    break;

	  case 's':
	    shadow_sw = 0;
	    break;

	  case 't':
	    if (trcdes) {
	      fprintf(stderr, "multiple trace count\n");
	      exit(1);
	    }
	    if (--argc == 0) {
	      fprintf(stderr, "missing specification of trace count\n");
	      exit(1);
	    }

	    {
	      char *descptr;

	      descptr = *++argv;
	      if ((trace_count = atouie(descptr, '\0', &descptr)) < 0) {
		fprintf(stderr, "invalid trace count %s\n", *argv);
		exit(1);
	      }
	    }
	    trcdes = TRUE;
	    break;

	  case 'v':
	    if (aspdes) {
	      fprintf(stderr, "multiple aspect ratio specification flag\n");
	      exit(1);
	    }
	    if (--argc == 0) {
	      fprintf(stderr, "missing specification of aspect ratio\n");
	      exit(1);
	    }

	    {
	      char *descptr = *++argv;

	      if ((aspect_h = atouie(descptr, ':', &descptr)) < 1) {
		fprintf(stderr, "invalid aspect_ratio description %s\n", *argv);
		exit(1);
	      }
	      if ((aspect_v = atouie(descptr, '\0', &descptr)) < 1) {
		fprintf(stderr, "invalid aspect_ratio description %s\n", *argv);
		exit(1);
	      }
	    }
	    aspdes = TRUE;
	    break;

	  default:
	    fprintf(stderr, "invalid flag %c\n", *flagptr);
	    exit(1);
	  }
      }
    }

    if (argc == 0) {
      fprintf(stderr, "missing source file description\n");
      exit(1);
    }
    strcpy(inpname, *argv);
    if(argc != 1) {
      fprintf(stderr, "%s may be not read\n", *++argv);
      exit(1);
    }

    if (!outdes) strcat(strcpy(outname, inpname), ".img");
    strcat(strcpy(maskout, inpname), ".msk");

    if (!aspdes) {
      aspect_h = reso_h;
      aspect_v = reso_v;
    }
  }

  if (super_level == 0) super_level = 1;

  shift_h = 0.0;
  shift_v = 0.0;
  screen_size_x = 35.0;
  screen_size_y = screen_size_x * aspect_v / aspect_h;
  pixel_size_x = screen_size_x / reso_h;
  pixel_size_y = screen_size_y / reso_v;
  screen_left_edge = screen_size_x * shift_h - screen_size_x / 2.0;
  screen_upper_edge = screen_size_y * shift_v + screen_size_y / 2.0;

  if ((fp = fopen(inpname, "r")) == NULL) {
    perror(inpname);
    exit(1);
  }

  clrpara();
  setpara(fp);
  fclose(fp);

  cam_refract();

  if ((wfp = fopen(outname, "w")) == NULL) {
    perror(outname);
    exit(1);
  }
  compute();
  fclose(wfp);
}

/* end */
