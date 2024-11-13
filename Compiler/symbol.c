#include <stdio.h>
#include <string.h>
#include "config.h"
#include "vector.h"
#include "object.h"
#include "symbol.h"
#include "subscan.h"

struct _rword {
  char rword_val;
  char *rwords;
};

static struct _rword rword[] = {
  SYMMOVE,		"move",
  SYMSCALE,		"scale",
  SYMROTATE,		"rotate",
  SYMCUBE,		"cube",
  SYMRECTANGLE,		"rectangle",
  SYMPLANE,		"plane",
  SYMSPHERE,		"sphere",
  SYMELLIPSOID,		"ellipsoid",
  SYMHYPERBOLOID,	"hyperboloid",
  SYMCONE,		"cone",
  SYMPARABOLOID,	"paraboloid",
  SYMCYLINDER,		"cylinder",
  SYMSUR2ND,		"sur2nd",
  SYMATTRIBUTE,		"attrib",
  SYMAMBIENCE,		"ambience",
  SYMLIGHT,		"light",
  SYMCAMERA,		"camera",
  SYMPUT,		"put",
  SYMDIFFUSIBILITY,	"diffusibility",
  SYMSPECULAR,		"specular",
  SYMTRANSPARENCY,	"transparency",
  SYMHIGHLIGHT,		"highlight",
  SYMREFRACT,		"refract",
  SYMPARALLEL,		"parallel",
  SYMINCLUDE,		"include",
  SYMTEXTURE,		"texture",
  SYMAT,		"at",
  SYMLUMINOSITY,	"luminosity",
  SYMPOSITION,		"position",
  SYMDIRECTION,		"direction",
  SYMTARGET,		"target",
  SYMZOOM,		"zoom",
  SYMPOINT,		"point",
  SYMFPOINT,		"fpoint",
  SYMSPOT,		"spot",
  SYMSKY,		"sky",
  SYMX,			"x",
  SYMY,			"y",
  SYMZ,			"z",
  NULL,			NULL
};

static SYMTAB *root;
static SYMTAB *symlast;
static int symdir;

static bool ismatch_symbol(symptr, s)
SYMTAB *symptr;
char *s;
{
  register int cmp;

  if ((cmp = strncmp(s, symptr->symword, MAXSYMLEN)) == 0) {
    symlast = symptr;
    return TRUE;
  }

  if (cmp < 0)
    if (symptr->lower_sym == NULL) {
      symlast = symptr;
      symdir = 0;
      return FALSE;
    }
    else return ismatch_symbol(symptr->lower_sym, s);
  else
    if (symptr->upper_sym == NULL) {
      symlast = symptr;
      symdir = 1;
      return FALSE;
    }
    else return ismatch_symbol(symptr->upper_sym, s);
}

void search_sym(s)
char *s;
{
  if (root)
    if (ismatch_symbol(root, s)) sym = symlast;
    else sym = NULL;
  else sym = NULL;
}

SYMTAB *def_symbol(s, type, value)
char *s;
int type;
union _symval *value;
{
  SYMTAB *symptr;
  int sl;

  if (root && ismatch_symbol(root, s)) return NULL;

  symptr = (SYMTAB *)malloc(sizeof(SYMTAB));
  if ((sl = strlen(s)) > MAXSYMLEN) sl = MAXSYMLEN;
  strncpy((symptr->symword = (char *)malloc(sl + 1)), s, sl);
  (symptr->symword)[sl] = '\0';

  switch (symptr->symtype = type)
    {
    case RESSYM:
      symptr->symval.code = value->code;
      break;

    case CSGSYM:
      symptr->symval.objptr = value->objptr;
      break;

    case ATTSYM:
      symptr->symval.atptr = value->atptr;
      break;

    default:
      int_error("illegal type symbol is defined");
    }

  symptr->lower_sym = symptr->upper_sym = NULL;

  if (root == NULL) root = symptr;
  else if (symdir == 0) {
    if (symlast->lower_sym != NULL) int_error("bad symbol link");
    else symlast->lower_sym = symptr;
  }
  else {
    if (symlast->upper_sym != NULL) int_error("bad symbol link");
    else symlast->upper_sym = symptr;
  }

  return symptr;
}

void init_sym()
{
   struct _rword *ptr;

   root = NULL;
   for (ptr = rword; ptr->rwords; ptr++) {
     union _symval temp;

     temp.code = ptr->rword_val;
     (void)def_symbol(ptr->rwords, RESSYM, &temp);
   }
}

/* end */
