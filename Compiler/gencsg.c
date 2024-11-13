#include <stdio.h>
#include "config.h"
#include "vector.h"
#include "object.h"
#include "matrix.h"
#include "tree.h"
#include "symbol.h"
#include "default.h"
#include "error.h"

#define LEFT 1
#define RIGHT 2

#define MAX_OPTFACT    20	/* max. number of factor of optibuf */
#define MAX_UNNOT_TERM 20	/* max. number of term in ~( ) */

typedef struct _optibuf {
   Cstruct *id;
   unsigned int openum;
   bool mposi[MAX_OPTFACT];
   bool mstatic[MAX_OPTFACT];
   char mmark[MAX_OPTFACT];
   SHAPE *member[MAX_OPTFACT];
   ATTR *mattrib[MAX_OPTFACT];
   struct _optibuf *next_opti;
} OPTI;

static OPTI *opti_root, *opti_last;
static OPTI *opti_new, *opti_back, *opti_next, *opti_hook;
static Cstruct *ctree_ptr;
static OPTI *unnot_list[MAX_UNNOT_TERM];
static int num_term;

static OPTI *copy_factor(dest, destfn, sorc, sorcfn)
OPTI *dest, *sorc;
int destfn, sorcfn;
{
  dest->mposi[destfn] = sorc->mposi[sorcfn];
  dest->mstatic[destfn] = sorc->mstatic[sorcfn];
  dest->member[destfn] = sorc->member[sorcfn];
  dest->mattrib[destfn] = sorc->mattrib[sorcfn];

  return dest;
}

static bool cmp_shape(arg1, arg2)
SHAPE *arg1, *arg2;
{
  if (arg1 == arg2) return TRUE;

  return
    arg1->func == arg2->func &&
    arg1->position.x == arg2->position.x &&
    arg1->position.y == arg2->position.y &&
    arg1->position.z == arg2->position.z &&
    arg1->shno == arg2->shno&&
    arg1->a == arg2->a &&
    arg1->b == arg2->b &&
    arg1->c == arg2->c &&
    arg1->d == arg2->d &&
    arg1->e == arg2->e &&
    arg1->f == arg2->f &&
    arg1->g == arg2->g &&
    arg1->h == arg2->h &&
    arg1->i == arg2->i &&
    arg1->j == arg2->j ;
}

static void optimize(idno)
Cstruct *idno;
{
  register OPTI *ptr, *refptr;

/*
  delete multiple factor ( a & a )
  delete term contains ( a & ~a )
*/
  for (ptr = opti_root; ptr; ) {
    opti_next = (opti_back = ptr)->next_opti;

    if (ptr->id == idno) {
      register unsigned int j;

      for (j = 0; j < ptr->openum; j++) {
	register unsigned int k;

	if (ptr->mattrib[j]->tranflag) continue;

	for (k = j + 1; k < ptr->openum; ) {
	  if (ptr->mattrib[k]->tranflag) {
	    k++;
	    continue;
	  }

	  if (cmp_shape(ptr->member[k], ptr->member[j]))
	    if (ptr->mposi[k] == ptr->mposi[j]) {
	      register int l;

	      for (l = k + 1; l < ptr->openum; l++)
		copy_factor(ptr, l - 1, ptr, l);
	      --(ptr->openum);
	    }
	    else {
	      ptr->openum = 0;
	      break;
	    }
	  else ++k;
	}
      }
    }

    if (ptr->openum == 0) {
      free((char *)ptr);

      if (opti_back == opti_last) opti_last = opti_hook;

      if (opti_back == opti_root) opti_root = opti_next;
      else opti_hook->next_opti = opti_next;
    }
    else opti_hook = opti_back;

    ptr = opti_next;
  }
/*
  delete term ( a & b & c ) if exists ( a & b ) etc.
*/
  for (refptr = opti_root; refptr; refptr = refptr->next_opti)
    if (refptr->id == idno) {
      register unsigned int refnum = refptr->openum;
      register unsigned int l;

      for (l = 0; l < refnum && !(refptr->mattrib[l]->tranflag); l++) ;
      if (l < refnum) continue;

      for (ptr = opti_root; ptr; ) {
	register unsigned int found_num = 0;

	opti_next = (opti_back = ptr)->next_opti;

	if (ptr->id == idno && ptr != refptr && ptr->openum >= refnum) {
	  while (found_num < refnum) {
	    register unsigned int l;
	    bool found = FALSE;

	    for (l = 0; l < ptr->openum && !(ptr->mattrib[l]->tranflag); l++)
	      if (cmp_shape(ptr->member[l], refptr->member[found_num]) &&
		  ptr->mposi[l] == refptr->mposi[found_num]) {
		found = TRUE;
		break;
	      }

	    if (found) ++found_num;
	    else break;
	  }
	}

	if (found_num == refnum) {
	  free((char *)ptr);
	  if (opti_back == opti_last) opti_last = opti_hook;
	  if (opti_back == opti_root) opti_root = opti_next;
	  else opti_hook->next_opti = opti_next;
	}
	else opti_hook = opti_back;
	ptr = opti_next;
      }
    }
}

static void connect_newopti()
{
  if (opti_root) opti_last->next_opti = opti_new;
  else opti_root = opti_new;

  opti_last = opti_new;
}

static OPTI *alloc_newopti()
{
  if ((opti_new = (OPTI *)malloc(sizeof(OPTI))) == NULL)
    error(1, TOO_OBJECT, NULL);

  opti_new->id = ctree_ptr;
  opti_new->openum = 0;
  opti_new->next_opti = NULL;

  return opti_new;
}

static OPTI *add_factor(ptr, type, shape)
OPTI *ptr;
int type;
char *shape;
{
  register unsigned int i;

  if (shape) {
    if ((i = ptr->openum++) >= MAX_OPTFACT) return 2;   /* ??? */
    ptr->mposi[i] = TRUE;
    ptr->mstatic[i] = (type != NEWSH_TYPE);
    ptr->member[i] = (SHAPE *)shape;
    ptr->mattrib[i] = dattrib;
  }

  return ptr;
}

static int make_or(side)
int side;
{
  int type;
  char *factor;

  if (side == LEFT) {   /* left operand */
    type = ctree_ptr->left_type;
    factor = ctree_ptr->left;
  }
  else {                /* right operand */
    type = ctree_ptr->right_type;
    factor = ctree_ptr->right;
  }

  switch (type)
    {
    case SHAPE_TYPE:
    case NEWSH_TYPE:
      add_factor(alloc_newopti(), type, factor);
      connect_newopti();
      break;

    case CSTRC_TYPE:
      {
	register OPTI *ptr;

	for (ptr = opti_root; ptr; ptr = ptr->next_opti)
	  if (ptr->id == (Cstruct *)factor) ptr->id = ctree_ptr;
      }
      break;

    default:
      return 1;
    }

  return 0;
}

static int unnots(termno)
int termno;
{
  if (termno < num_term) {
    OPTI *where;
    int who, num_fact;

    num_fact = (where = unnot_list[termno])->openum;
    for (who = 0; who < num_fact; who++) {
      if (who > 0) {
	register int i;

	alloc_newopti();
	for (i = 0; i < termno; i++)
	  copy_factor(opti_new, i, opti_last, i);
      }
      copy_factor(opti_new, termno, where, who);
      opti_new->mposi[termno] = !(opti_new->mposi[termno]);
      if (unnots(termno + 1)) return 1;
    }
  }
  else {
    opti_new->openum = num_term;
    connect_newopti();
  }

  return 0;
}

static int unnot()
{
  alloc_newopti();
  return unnots(0);
}

static OPTI *ins_factor(ptr, type, shape)
OPTI *ptr;
char type;
char *shape;
{
  register unsigned int i, j, k;

  if (shape) {
    if ((i = ptr->openum++) >= MAX_OPTFACT) return 2;   /* ??? */
    for (j = i; j; j--) {
      k = j - 1;
      ptr->mposi[j] = ptr->mposi[k];
      ptr->mstatic[j] = ptr->mstatic[k];
      ptr->member[j] = ptr->member[k];
      ptr->mattrib[j] = ptr->mattrib[k];
    }
    ptr->mposi[0] = TRUE;
    ptr->mstatic[0] = (type != NEWSH_TYPE);
    ptr->member[0] = (SHAPE *)shape;
    ptr->mattrib[0] = dattrib;
  }

  return ptr;
}

static OPTI *cat_term(dest, sorc)
OPTI *dest, *sorc;
{
  register int j;

  for (j = 0; j < sorc->openum; j++) {
    register int i;

    if ((i = dest->openum++) >= MAX_OPTFACT) return 2;   /* ??? */
    copy_factor(dest, i, sorc, j);
  }

  return dest;
}

static void clr_mark(id)
Cstruct *id;
{
  register OPTI *ptr;
  register unsigned int i;

  for (ptr = opti_root; ptr; ptr = ptr->next_opti)
    if (ptr->id == id)
      for (i = 0; i < ptr->openum; i++)
	ptr->mmark[i] = 0;
}

static void chg_shptr(id, old_shptr, new_shptr)
Cstruct *id;
SHAPE *old_shptr, *new_shptr;
{
  register OPTI *ptr;
  register unsigned int i;

  for (ptr = opti_root; ptr; ptr = ptr->next_opti)
    if (ptr->id == id)
      for (i = 0; i < ptr->openum; i++)
	if (ptr->member[i] == old_shptr) {
	  ptr->member[i] = new_shptr;
	  ptr->mmark[i] = 1;
	}
}

static SHAPE *rot_surface(sh, m)
SHAPE *sh;
matrixp m;
{
  double a, b, c, d, e, f, g, h, i;

  a = sh->a;
  b = sh->b;
  c = sh->c;
  d = sh->d;
  e = sh->e;
  f = sh->f;
  g = sh->g;
  h = sh->h;
  i = sh->i;

  sh->a = a * m[0][0] * m[0][0]
        + b * m[1][0] * m[1][0]
	+ c * m[2][0] * m[2][0]
	+ d * m[1][0] * m[2][0]
	+ e * m[0][0] * m[2][0]
	+ f * m[0][0] * m[1][0];

  sh->b = a * m[0][1] * m[0][1]
        + b * m[1][1] * m[1][1]
	+ c * m[2][1] * m[2][1]
	+ d * m[1][1] * m[2][1]
	+ e * m[0][1] * m[2][1]
	+ f * m[0][1] * m[1][1];

  sh->c = a * m[0][2] * m[0][2]
        + b * m[1][2] * m[1][2]
	+ c * m[2][2] * m[2][2]
	+ d * m[1][2] * m[2][2]
	+ e * m[0][2] * m[2][2]
	+ f * m[0][2] * m[1][2];

  sh->d = 2.0 * (a * m[0][1] * m[0][2]
	       + b * m[1][1] * m[1][2]
	       + c * m[2][1] * m[2][2])
        + d * (m[1][1] * m[2][2] + m[2][1] * m[1][2])
	+ e * (m[0][1] * m[2][2] + m[2][1] * m[0][2])
	+ f * (m[0][1] * m[1][2] + m[1][1] * m[0][2]);

  sh->e = 2.0 * (a * m[0][0] * m[0][2]
	       + b * m[1][0] * m[1][2]
	       + c * m[2][0] * m[2][2])
        + d * (m[1][0] * m[2][2] + m[2][0] * m[1][2])
	+ e * (m[0][0] * m[2][2] + m[2][0] * m[0][2])
	+ f * (m[0][0] * m[1][2] + m[1][0] * m[0][2]);

  sh->f = 2.0 * (a * m[0][0] * m[0][1]
	       + b * m[1][0] * m[1][1]
	       + c * m[2][0] * m[2][1])
        + d * (m[1][0] * m[2][1] + m[2][0] * m[1][1])
	+ e * (m[0][0] * m[2][1] + m[2][0] * m[0][1])
	+ f * (m[0][0] * m[1][1] + m[1][0] * m[0][1]);

  sh->g = g * m[0][0] + h * m[1][0] + i * m[2][0];
  sh->h = g * m[0][1] + h * m[1][1] + i * m[2][1];
  sh->i = g * m[0][2] + h * m[1][2] + i * m[2][2];

  return sh;
}

SHAPE *nplane(sh)
SHAPE *sh;
{
  vector v;
  extern bool normalize();

  v.x = sh->g;
  v.y = sh->h;
  v.z = sh->i;
  normalize(&v);
  sh->g = v.x;
  sh->h = v.y;
  sh->i = v.z;

  return sh;
}

static SHAPE *sca_surface(sh, m)
SHAPE *sh;
matrixp m;
{
  sh->a /= m[0][0] * m[0][0];
  sh->b /= m[1][1] * m[1][1];
  sh->c /= m[2][2] * m[2][2];
  sh->d /= m[1][1] * m[2][2];
  sh->e /= m[0][0] * m[2][2];
  sh->f /= m[0][0] * m[1][1];
  sh->g /= m[0][0];
  sh->h /= m[1][1];
  sh->i /= m[2][2];

  if (sh->func == PLANE) nplane(sh);

  return sh;
}

int genCSG()
{
  register int i;

  ctree_ptr = ctree;
  opti_root = NULL;

  for (i = 0; i < ctree_count; i++)
    {
      switch (ctree_ptr->operator)
	{
	case OP_EQUAL:
	  if (ctree_ptr->left_type != LVALUE_TYPE) int_error("not lvalue");
	  if (make_or(RIGHT)) return 1;
	  {
	    ORLIST *orptr, *or_root, *or_last;
	    OBJECT *object;
	    register OPTI *ptr;

	    or_root = or_last = NULL;

	    for (ptr = opti_root; ptr; ptr = ptr->next_opti)
	      if (ptr->id == ctree_ptr && ptr->openum > 0) {
		register unsigned int k, num_factor;
		ANDFACT *part, *last_andfact;

		if ((orptr = (ORLIST *)malloc(sizeof(ORLIST))) == NULL)
		  error(1, TOO_OBJECT, NULL);
		orptr->next_term = NULL;

		if (or_last) or_last->next_term = orptr;
		else or_root = orptr;

		num_factor = ptr->openum;
		if ((orptr->and_root = part = (ANDFACT *)calloc(num_factor, sizeof(ANDFACT))) == NULL)
		  error(1, TOO_OBJECT, NULL);

		for (k = 0; k < num_factor; k++) {
		  part->shape = ptr->member[k];
		  part->positive = ptr->mposi[k];
		  part->attribute = ptr->mattrib[k];
		  ptr->mstatic[k] = TRUE;
		  last_andfact = part;
		  last_andfact->next_factor = ++part;
		}
		last_andfact->next_factor = NULL;
		or_last = orptr;
	      }

	    if ((object = (OBJECT *)malloc(sizeof(OBJECT))) == NULL)
	      error(1, TOO_OBJECT, NULL);

	    object->or_root = or_root;
	    ((SYMTAB *)(ctree_ptr->left))->symval.objptr = object;
	  }
	  break;

	case OP_UNION:
	  if (make_or(LEFT)) return 1;
	  if (make_or(RIGHT)) return 1;
	  optimize(ctree_ptr);
	  break;

	case OP_INTER:
	  {
	    int mode = 0;

	    switch (ctree_ptr->left_type)
	      {
	      case CSTRC_TYPE:
		mode |= 2;
		break;

	      case SHAPE_TYPE:
	      case NEWSH_TYPE:
		break;

	      default:
		return 1;
	      }

	    switch (ctree_ptr->right_type)
	      {
	      case CSTRC_TYPE:
		mode |= 1;
		break;

	      case SHAPE_TYPE:
	      case NEWSH_TYPE:
		break;

	      default:
		return 1;
	      }

	    switch (mode)
	      {
	      case 0:   /* ( operand & operand ) */
		add_factor(add_factor(alloc_newopti(), ctree_ptr->left_type, ctree_ptr->left), ctree_ptr->right_type, ctree_ptr->right);
		connect_newopti();
		break;

	      case 1:   /* ( operand & C-tree ) */
	      case 2:   /* ( C-tree & operand ) */
		{
		  Cstruct *dest = (Cstruct *)
		    (mode == 1 ? ctree_ptr->right : ctree_ptr->left);
		  OPTI *ptr;

		  for (ptr = opti_root; ptr; ptr = ptr->next_opti)
		    if (ptr->id == dest) {
		      if (mode == 1)
			ins_factor(ptr, ctree_ptr->left_type, ctree_ptr->left);
		      else
			add_factor(ptr, ctree_ptr->right_type, ctree_ptr->right);
		      ptr->id = ctree_ptr;
		    }
		}
		break;

	      case 3:   /* ( C-tree & C-tree ) */
		{
		  register OPTI *ptrL, *ptrR;
		  Cstruct *what_Lnode, *what_Rnode;

		  what_Lnode = (Cstruct *)(ctree_ptr->left);
		  what_Rnode = (Cstruct *)(ctree_ptr->right);

		  for (ptrL = opti_root; ptrL && ptrL->id != ctree_ptr; ) {
		    opti_next = (opti_back = ptrL)->next_opti;

		    if (ptrL->id == what_Lnode) {

		      for (ptrR = opti_root;
			   ptrR && ptrR->id != ctree_ptr;
			   ptrR = ptrR->next_opti)
			if (ptrR->id == what_Rnode) {
			  cat_term(cat_term(alloc_newopti(), ptrL), ptrR);
			  connect_newopti();
			}

		      free((char *)ptrL);

		      if (opti_back == opti_last) opti_last = opti_hook;

		      if (opti_back == opti_root) opti_root = opti_next;
		      else opti_hook->next_opti = opti_next;
		    }
		    else opti_hook = opti_back;

		    ptrL = opti_next;
		  }

		  for (ptrR = opti_root; ptrR && ptrR->id != ctree_ptr; ) {
		    opti_next = (opti_back = ptrR)->next_opti;
		    if (ptrR->id == what_Rnode) {
		      free((char *)ptrR);

		      if (opti_back == opti_last) opti_last = opti_hook;

		      if (opti_back == opti_root) opti_root = opti_next;
		      else opti_hook->next_opti = opti_next;
		    }
		    else opti_hook = opti_back;

		    ptrR = opti_next;
		  }
		}
		break;
	      }
	  }
	  optimize(ctree_ptr);
	  break;

	case OP_NOT:
	  switch (ctree_ptr->right_type)
	    {
	    case SHAPE_TYPE:
	    case NEWSH_TYPE:
	      add_factor(alloc_newopti(), ctree_ptr->right_type, ctree_ptr->right);
	      opti_new->mposi[0] = FALSE;
	      connect_newopti();
	      break;

	    case CSTRC_TYPE:
	      {
		register OPTI *ptr;
		Cstruct *what_node;
		OPTI **unlist_ptr;

		what_node = (Cstruct *)(ctree_ptr->right);
		unlist_ptr = unnot_list;
		num_term = 0;

		for (ptr = opti_root; ptr; ptr = ptr->next_opti) {
		  if (ptr->id == what_node) {
		    if (num_term >= MAX_UNNOT_TERM) return 1;
		    *unlist_ptr++ = ptr;
		    ++num_term;
		  }
		}

		if (num_term == 0)
		  int_error("undefined C-structure as operand of NOT operation");
		if (unnot()) return 1;

		for (ptr = opti_root; ptr && ptr->id != ctree_ptr; ) {
		  opti_next = (opti_back = ptr)->next_opti;
		  if (ptr->id == what_node) {
		    free((char *)ptr);

		    if (opti_back == opti_last) opti_last = opti_hook;

		    if (opti_back == opti_root) opti_root = opti_next;
		    else opti_hook->next_opti = opti_next;
		  }
		  else opti_hook = opti_back;

		  ptr = opti_next;
		}
	      }
	      break;

	    default:
	      return 1;
	    }
	  optimize(ctree_ptr);
	  break;

	case ADDATR:
	  {
	    register OPTI *ptr;
	    register ATTR *at;

	    if (ctree_ptr->left_type != ATTR_TYPE) error(1, BAD_SYNTAX, NULL);

	    at = (ATTR *)ctree_ptr->left;

	    if (make_or(RIGHT)) error(1, BAD_SYNTAX, NULL);

	    for (ptr = opti_root; ptr; ptr = ptr->next_opti)
	      if (ptr->id == ctree_ptr) {
		register unsigned int l, n;

		n = ptr->openum;
		for (l = 0; l < n; l++)
		  if (ptr->mattrib[l] == dattrib) ptr->mattrib[l] = at;
	      }
	  }
	  break;

	case SYMMOVE:
	case SYMSCALE:
	case SYMROTATE:
	  {
	    register OPTI *ptr;
	    matrixp m;

	    if (ctree_ptr->right_type != MATRIX_TYPE) return 1;

	    if (make_or(LEFT)) error(1, BAD_SYNTAX, NULL);

	    clr_mark(ctree_ptr);
	    m = (matrixp)ctree_ptr->right;
	    for (ptr = opti_root; ptr; ptr = ptr->next_opti) {
	      if (ptr->id == ctree_ptr) {
		register unsigned int l;

		for (l = 0; l < ptr->openum; l++) {
		  if (ptr->mmark[l] == 0) {
		    SHAPE *new_shape;
		    extern SHAPE *alloc_shape(), *copy_shape();

		    if (ptr->mstatic[l]) {
		      new_shape = copy_shape(alloc_shape(), ptr->member[l]);
		      new_shape->shno = 0;
		    }
		    else new_shape = ptr->member[l];

		    switch (ctree_ptr->operator)
		      {
		      case SYMSCALE:
			sca_surface(new_shape, m);
			break;

		      case SYMROTATE:
			rot_surface(new_shape, m);
			break;
		      }
		    apply(&(new_shape->position), m);
		    chg_shptr(ctree_ptr, ptr->member[l], new_shape);
		  }
		}
	      }
	    }
	    free((char *)m);
	  }
	  break;

	default:
	  return 3;   /* bad operator */
	}

      ++ctree_ptr;
    }

  return 0;
}

void free_opti()
{
  while (opti_root) {
    opti_next = opti_root->next_opti;
    free((char *)opti_root);
    opti_root = opti_next;
  }
}

/* end */
