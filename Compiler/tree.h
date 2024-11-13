typedef struct _C_struct {
  char operator;
  char left_type;
  char right_type;
  char *left;
  char *right;
} Cstruct;

extern Cstruct *ctree;
extern int ctree_count;

#define	NULL_TYPE	0
#define	CSTRC_TYPE	1
#define	LVALUE_TYPE	2
#define	NEWSH_TYPE	3
#define	SHAPE_TYPE	4
#define	ATTR_TYPE	5
#define	MATRIX_TYPE	6

extern void alloc_ctree();
extern void init_expr();
extern int ent_op();
extern int ent_CSGshape();
extern int ent_operand();
extern int genCSG();
extern void free_opti();
