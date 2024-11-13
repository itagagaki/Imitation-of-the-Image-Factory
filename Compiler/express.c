#include <stdio.h>
#include "config.h"
#include "vector.h"
#include "object.h"
#include "tree.h"
#include "symbol.h"
#include "default.h"

#define MAXOPS	       20	/* max. operator srack level */
#define MAXCTREE      100	/* max. number of C-structure tree ( nearly equal number of operator ) */

#define Un 1	/* cut () */
#define Lo 2	/* push */
#define Up 3	/* tree out last terminal */
#define Er 4	/* error - bad syntax */
#define Uc 5	/* error - found ; expcted close of ( */
#define Si 6
#define Di 7	/* replace - to &~ */

static char prec_table[16][16] = {
/*	   BOE =  += &= -= +  &  ~  -  (  )  ;  AT mv sc rt */
/* BOE */ { Er,Lo,Lo,Lo,Lo,Lo,Lo,Lo,Di,Lo,Er,Si,Lo,Lo,Lo,Lo },
/*  =  */ { Er,Lo,Lo,Lo,Lo,Lo,Lo,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  += */ { Er,Lo,Lo,Lo,Lo,Lo,Lo,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  &= */ { Er,Lo,Lo,Lo,Lo,Lo,Lo,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  -= */ { Er,Lo,Lo,Lo,Lo,Lo,Lo,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  +  */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  &  */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  ~  */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  -  */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Lo,Lo,Lo,Lo },
/*  (  */ { Er,Lo,Lo,Lo,Lo,Lo,Lo,Lo,Di,Lo,Un,Uc,Lo,Lo,Lo,Lo },
/*  )  */ { Er,Up,Up,Up,Up,Up,Up,Up,Di,Er,Up,Up,Up,Up,Up,Up },
/*  ;  */ { Er,Er,Er,Er,Er,Er,Er,Er,Di,Er,Er,Er,Er,Er,Er,Er },
/* ATR */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Er,Up,Up,Up },
/* mov */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Er,Up,Up,Up },
/* sca */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Er,Up,Up,Up },
/* rot */ { Er,Up,Up,Up,Up,Up,Up,Lo,Di,Lo,Up,Up,Er,Up,Up,Up }
};

struct _opstack {
   char operator;
   char ope_type;
   char *operand;	/* pointer for SHAPE or Cstruct */
} opstack[MAXOPS];

Cstruct *ctree, *ctree_ptr;

int ops_ptr;
SYMTAB *pending_symbol;
int ctree_num, ctree_free, ctree_count;

int ent_op();

void alloc_ctree()
{
  if ((ctree = (Cstruct *)calloc( MAXCTREE, sizeof(Cstruct))) == NULL) {
    fprintf(stderr, "No enough memory space for C-structure tree\n");
    exit(1);
  }

  ctree_num = MAXCTREE;
}

void init_expr()
{
  ctree_ptr = ctree;
  ctree_free = ctree_num;
  ctree_count = 0;
  pending_symbol = (SYMTAB *)NULL;
  ops_ptr = 0;
  if (ent_op(BOE)) int_error("unable initialize expression analyzer");
}

int ent_CSGshape(symptr)
SYMTAB *symptr;
{
  if (pending_symbol) return 1;
  pending_symbol = symptr;
  return 0;
}

static int exp_op(op)
int op;
{
  int finished = 0;

  while (ops_ptr)
    {
      switch(prec_table[opstack[ops_ptr - 1].operator][op])
	{
	case Lo:
	  finished = 1;
	  break;

	case Un:
	  --ops_ptr;
	  if (opstack[ops_ptr - 1].ope_type != NULL_TYPE) return 1;  /* bad syntax */ 
	  opstack[ops_ptr - 1].ope_type = opstack[ops_ptr].ope_type;
	  opstack[ops_ptr - 1].operand = opstack[ops_ptr].operand;
	  opstack[ops_ptr].operator = NOOP;
	  finished = 2;
	  break;

	case Up:
	  {
	    int top;
	    char left_type, right_type;
	    char *left, *right;

	    --ops_ptr;
	    top = opstack[ops_ptr].operator;
	    left_type = opstack[ops_ptr - 1].ope_type;
	    left = opstack[ops_ptr - 1].operand;
	    right_type = opstack[ops_ptr].ope_type;
	    right = opstack[ops_ptr].operand;

	    if (top == NOOP) return 1;

	    if (top == OP_EQUNION || top == OP_EQINTER || top == OP_EQDIFF)
	      {
		ent_op(OP_EQUAL);
		ent_CSGshape((SYMTAB *)left);

		switch (top)
		  {
		  case OP_EQUNION:
		    ent_op(OP_UNION);
		    break;

		  case OP_EQINTER:
		    ent_op(OP_INTER);
		    break;

		  case OP_EQDIFF:
		    ent_op(OP_DIFF);
		    break;
		  }

		opstack[--ops_ptr].ope_type = right_type;
		opstack[ops_ptr++].operand = right;
		ent_op(op);
		finished = 2;
		break;
	      }

	    if (top == OP_EQUAL) {
	      if (left_type != LVALUE_TYPE) return 5;  /* illegal = operator */
	      if (right_type == LVALUE_TYPE) return 7;
	    }
	    else {
	      if (left_type == LVALUE_TYPE) return 7;
	      if (right_type == LVALUE_TYPE) return 7;
	    }

	    if (left_type == NULL_TYPE && top != OP_NOT)
	      return 1;   /* bad syntax (no left operand in not ~ operation) */

	    if (right_type == NULL_TYPE)
	      return 1;   /* bad syntax (no right operand) */

	    if (ctree_free == 0)
	      return 3;   /* too many C-structure */

	    ctree_ptr->operator = top;
	    ctree_ptr->left_type = left_type;
	    ctree_ptr->left = left;
	    ctree_ptr->right_type = right_type;
	    ctree_ptr->right = right;
	    opstack[--ops_ptr].ope_type = CSTRC_TYPE;
	    opstack[ops_ptr++].operand = (char *)ctree_ptr++;
	    ++ctree_count;
	    --ctree_free;
	  }
	  break;

	case Si:
	  if (opstack[ops_ptr - 1].ope_type != CSTRC_TYPE) return 6;  /* bad syntax */
	  finished = 2;
	  break;

	case Di:
	  if (ent_op(OP_INTER)) return 1;
	  if (ent_op(OP_NOT)) return 1;
	  finished = 2;
	  break;

	case Uc:
	  return 4;   /* found ; expected close of( */

	case Er:
	default:
	  return 1;   /* bad syntax */
	}

      if (finished) break;
    }

  if (finished != 2) {
    if(ops_ptr >= MAXOPS) return 2;   /* overflow of operator stack */
    opstack[ops_ptr].operator = op;
    opstack[ops_ptr++].ope_type = NULL_TYPE;
  }

  return 0;
}

int ent_operand(type, operand)
int type;
char *operand;
{
  register int opop;

  if (pending_symbol) return 1;

  opop = ops_ptr - 1;
  if (opstack[opop].operator == NOOP || opstack[opop].ope_type != NULL_TYPE)
    return 1;   /* bad syntax */

  opstack[opop].ope_type = type;
  opstack[opop].operand = operand;

  return 0;
}

static int exp_object(symptr)
SYMTAB *symptr;
{
  ORLIST *object;
  bool isroot_object = TRUE;

  if (object = symptr->symval.objptr->or_root)
    {
      ent_op(SYMLPAR);

      do {
	ANDFACT *part = object->and_root;
	bool isroot_part = TRUE;

	if (isroot_object) isroot_object = FALSE;
	else {
	  ent_op(OP_UNION);
	  object = object->next_term;
	}

	do {
	  if (isroot_part) isroot_part = FALSE;
	  else {
	    ent_op(OP_INTER);
	    part = part->next_factor;
	  }

	  if (!(part->positive)) ent_op(OP_NOT);

	  if (part->attribute != dattrib) {
	    ent_operand(ATTR_TYPE, (char *)( part->attribute));
	    ent_op(ADDATR);
	  }

	  ent_operand(SHAPE_TYPE, (char *)( part->shape));

	} while(part->next_factor);

      } while(object->next_term);

      ent_op(SYMRPAR);
    }
  else ent_operand(SHAPE_TYPE, (char *)NULL);

  return 0;
}

/*
enter operator to stack and analize with operator precedence grammer

return:
	0	no error
	1	bad syntax
	2	overflow of operator stack
	3	too many C-structure
	4	found ; expected close of (
	5	illegal = operator
*/
int ent_op(op)
int op;
{
  SYMTAB *symptr;

  if (pending_symbol) {
    symptr = pending_symbol;
    pending_symbol = (SYMTAB *)NULL;

    if (op == OP_EQUAL || op == OP_EQUNION || op == OP_EQINTER || op == OP_EQDIFF)
      ent_operand(LVALUE_TYPE, (char *)symptr);
    else
      exp_object(symptr);
   }

  return exp_op(op);
}

#if 0

bool isequop(op)
int op;
{
  return (op == OP_EQUAL
	  || op == OP_EQUNION
	  || op == OP_EQINTER
	  || op == OP_EQDIFF
	  );
}
#endif

#if 0

char *opchr(c)
int c;
{
  static char *ss[] = {
    "  [   ",
    "  =   ",
    "  +=  ",
    "  &=  ",
    "  -=  ",
    "  +   ",
    "  &   ",
    "  ~   ",
    "  -   ",
    " (  ",
    "  )  ",
    " move ",
    "rotate",
    "scale ",
    " atsh ",
    "  ;   "
    };

  return c < 16 ? ss[c] : "  .   ";
}

void print_tree()
{
  Cstruct *ctree_ptr;
  register int i;

  ctree_ptr = ctree;

  for(i = 0; i < ctree_count; i++)
    {
      printf("[%x] : ", ctree_ptr);
      switch(ctree_ptr->left_type)
	{
	case CSTRC_TYPE:
	  printf("[%x]", ctree_ptr->left);
	  break;

	case SHAPE_TYPE:
	case NEWSH_TYPE:
	  printf("%x", ctree_ptr->left);
	  break;

	case ATTR_TYPE:
	  printf("ATTR %x", ctree_ptr->left);
	  break;

	case LVALUE_TYPE:
	  printf("(%x)", ctree_ptr->left);
	  break;

	default:
	  printf("NULL");
	  break;
	}

      printf(" %s ", opchr(ctree_ptr->operator));

      switch(ctree_ptr->right_type)
	{
	case CSTRC_TYPE:
	  printf("[%x]", ctree_ptr->right);
	  break;

	case SHAPE_TYPE:
	case NEWSH_TYPE:
	  printf("%x", ctree_ptr->right);
	  break;

	case ATTR_TYPE:
	  printf("ATTR %x", ctree_ptr->right);
	  break;

	case LVALUE_TYPE:
	  printf("(%x)", ctree_ptr->right);
	  break;

	default:
	  printf("NULL");
	  break;
	}

      putchar('\n');
      ++ctree_ptr;
    }
}

void print_stack()
{
  int i;

  putchar('\n');

  for (i = ops_ptr - 1; i >= 0; i--)
    {
      printf("%d : ", i);

      switch(opstack[i].operator)
	{
	case NOOP:
	  printf("  NOOP  ");
	  break;

	default:
	  printf(" %s ", opchr(opstack[i].operator));
	  break;
	}

      printf(" | ");

      switch(opstack[i].ope_type)
	{
	case NULL_TYPE:
	  puts("NULL");
	  break;

	case SHAPE_TYPE:
	  printf("SHAPE  %x\n", opstack[i].operand);
	  break;

	case NEWSH_TYPE:
	  printf("NEWSH  %x\n", opstack[i].operand);
	  break;

	case ATTR_TYPE:
	  printf("ATTR   %x\n", opstack[i].operand);
	  break;

	case MATRIX_TYPE:
	  printf("MATRIX %x\n", opstack[i].operand);
	  break;

	case CSTRC_TYPE:
	  printf("CTRRC  %x\n", opstack[i].operand);
	  break;

	case LVALUE_TYPE:
	  printf("LVALUE %x\n", opstack[i].operand);
	  break;

	defalut:
	  puts("??????");
	  break;
	}
    }
}

#endif

/* end */
