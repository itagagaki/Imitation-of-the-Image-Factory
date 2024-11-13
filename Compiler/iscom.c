/*
	ISCOM image score compiler 1.0
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "config.h"
#include "form.h"
#include "math2.h"
#include "vector.h"
#include "object.h"
#include "light.h"
#include "intcode.h"
#include "matrix.h"
#include "tree.h"
#include "symbol.h"
#include "token.h"
#include "subscan.h"
#include "default.h"
#include "error.h"
#include "output.h"

#define FNAM_LEN 80

FILE *wfp;		/* creating FILE pointer */
char *outfile;		/* pointer for name of creating file */
char *inpfile;		/* pointer for name of current input file */
int err_count;

static int shno, atno;
int camno, lino;
bool sky_def;
ATTR *dattrib;


static ATTR *alloc_attr()
{
  ATTR *a;

  if ((a = (ATTR *)malloc(sizeof(ATTR))) == NULL) error(1, TOO_ATTR, NULL);
  a->atno = 0;
  return a;
}

static ATTR *copy_attr(dest, src)
ATTR *dest, *src;
{
  dest->atno = src->atno;
  dest->tranflag = src->tranflag;
  dest->diff_red = src->diff_red;
  dest->diff_grn = src->diff_grn;
  dest->diff_blu = src->diff_blu;
  dest->spec_red = src->spec_red;
  dest->spec_grn = src->spec_grn;
  dest->spec_blu = src->spec_blu;
  dest->tran_red = src->tran_red;
  dest->tran_grn = src->tran_grn;
  dest->tran_blu = src->tran_blu;
  dest->lumi_red = src->lumi_red;
  dest->lumi_grn = src->lumi_grn;
  dest->lumi_blu = src->lumi_blu;
  dest->refract = src->refract;
  dest->hilit = src->hilit;
  dest->hilit_diff = src->hilit_diff;
  dest->texmap = src->texmap;
  return dest;
}

static ATTR *fget_attribute()
{
  ATTR *at;
  union _symval temp;
  unsigned char token;
  double str, red, grn, blu;

  if (lexi() != SYMTOK) error(1, BAD_SYNTAX, NULL);

  temp.atptr = at = copy_attr(alloc_attr(), dattrib);
  def_symbol(word, ATTSYM, &temp);

  check_token(OP_EQUAL);
  check_token(SYMLMPAR);

  for (;;) {
    lexi();

    if (toktype == OP_TOK && opcode == SYMRMPAR) break;

    if (toktype != SYMTOK || sym == NULL || sym->symtype != RESSYM)
      error(1, BAD_SYNTAX, NULL);

    token = sym->symval.code;
    check_token(OP_EQUAL);

    switch (token) {

    case SYMDIFFUSIBILITY:
      finput(&str, SYMLPAR, TRUE);
      finput(&red, SYMCOMMA, TRUE);
      finput(&grn, SYMCOMMA, TRUE);
      finput(&blu, SYMRPAR, TRUE);
      check_token(SYMEXPTERM);
      at->diff_red = red * str;
      at->diff_grn = grn * str;
      at->diff_blu = blu * str;
      break;

    case SYMSPECULAR:
      finput(&str, SYMLPAR, TRUE);
      finput(&red, SYMCOMMA, TRUE);
      finput(&grn, SYMCOMMA, TRUE);
      finput(&blu, SYMRPAR, TRUE);
      check_token(SYMEXPTERM);
      at->spec_red = red * str;
      at->spec_grn = grn * str;
      at->spec_blu = blu * str;
      break;

    case SYMHIGHLIGHT:
      finput(&(at->hilit), SYMCOMMA, TRUE);
      finput(&(at->hilit_diff), SYMEXPTERM, TRUE);
      break;

    case SYMTRANSPARENCY:
      finput(&str, SYMLPAR, TRUE);
      finput(&red, SYMCOMMA, TRUE);
      finput(&grn, SYMCOMMA, TRUE);
      finput(&blu, SYMRPAR, TRUE);
      check_token(SYMEXPTERM);
      at->tran_red = 1.0 - red * str;
      at->tran_grn = 1.0 - grn * str;
      at->tran_blu = 1.0 - blu * str;
      if (at->tran_red < 1.0 || at->tran_grn < 1.0 || at->tran_blu < 1.0)
	at->tranflag = TRUE;
      break;

    case SYMREFRACT:
      finput(&(at->refract), SYMEXPTERM, TRUE);
      break;

    case SYMLUMINOSITY:
      finput(&str, SYMLPAR, TRUE);
      finput(&red, SYMCOMMA, TRUE);
      finput(&grn, SYMCOMMA, TRUE);
      finput(&blu, SYMRPAR, TRUE);
      check_token(SYMEXPTERM);
      at->lumi_red = red * str;
      at->lumi_grn = grn * str;
      at->lumi_blu = blu * str;
      break;

    case SYMTEXTURE:
      {
	double tmp;

	finput(&tmp, SYMEXPTERM, TRUE);
	at->texmap = tmp;
      }
      break;

    default:
      error(1, BAD_SYNTAX, NULL);
    }
  }

  return at;
}

static vector *get_vector(v, omissible)
vector *v;
bool omissible;
{
  check_token(SYMLBRACKET);
  finput(&(v->x), SYMCOMMA,    omissible);
  finput(&(v->y), SYMCOMMA,    omissible);
  finput(&(v->z), SYMRBRACKET, omissible);
  return v;
}

static matrixp get_move()
{
   vector v;

   v.x = v.y = v.z = 0.0;
   get_vector(&v, TRUE);
   check_token(SYMRPAR);
   return move(unit(alloc_matrix()), &v);
}

static matrixp get_scale()
{
   vector v;

   v.x = v.y = v.z = 1.0;
   get_vector(&v, TRUE);
   check_token(SYMRPAR);
   return scale(unit(alloc_matrix()), &v);
}

static matrixp get_rotate()
{
  matrixp m;
  vector origin, axis, mmv;
  double angle, l;

  get_vector(&origin, FALSE);
  switch (fget_token()) {

  case SYMAXIS:
    get_vector(&axis, FALSE);
    axis.x -= origin.x;
    axis.y -= origin.y;
    axis.z -= origin.z;
    break;

  case SYMOAXIS:
    axis.x = axis.y = axis.z = 0.0;
    get_vector(&axis, TRUE);
    break;

  default:
    error(1, BAD_SYNTAX, NULL);
  }

  check_token(SYMCOMMA);
  finput(&angle, SYMRPAR, FALSE);
  angle = DEGtoRAD(angle);

  if (( l = vlen(&axis)) == 0.0) error(1, BAD_AXIS, NULL);

  m = unit(alloc_matrix());
  mmv.x = -(origin.x);
  mmv.y = -(origin.y);
  mmv.z = -(origin.z);
  move(m, &mmv);

  if      (axis.y == 0.0 && axis.z == 0.0) rotx(m, sign(axis.x) * angle);
  else if (axis.x == 0.0 && axis.z == 0.0) roty(m, sign(axis.y) * angle);
  else if (axis.x == 0.0 && axis.y == 0.0) rotz(m, sign(axis.z) * angle);
  else {
    double si, ci, sj, cj, v;

    v = sqrt(axis.y * axis.y + axis.z * axis.z);
    si = axis.y / v;
    ci = axis.z / v;
    sj = axis.x / l;
    cj = v / l;
    rotxp(m, si, ci);
    rotyp(m, -sj, cj);
    rotz(m, angle);
    rotyp(m, sj, cj);
    rotxp(m, -si, ci);
   }

   move(m, &origin);

   return m;
}

static void affin(func)
int func;
{
  matrixp m;

  check_token(SYMLPAR);
  ent_op(func);

  switch(func) {
  case SYMMOVE:

    m = get_move();
    break;

  case SYMSCALE:
    m = get_scale();
    break;

  case SYMROTATE:
    m = get_rotate();
    break;
  }

  ent_operand(MATRIX_TYPE, (char *)m);
}

SHAPE *alloc_shape()
{
  SHAPE *sh;

  if ((sh = (SHAPE *)malloc(sizeof(SHAPE))) == NULL) error(1, TOO_SHAPE, NULL);

  sh->shno = 0;

  sh->position.x =
  sh->position.y =
  sh->position.z =
  sh->a =
  sh->b =
  sh->c =
  sh->d =
  sh->e =
  sh->f =
  sh->g =
  sh->h =
  sh->i =
  sh->j = 0.0;
  return sh;
}

SHAPE *copy_shape(dest, srce)
SHAPE *dest, *srce;
{
   dest->shno = srce->shno;
   dest->func = srce->func;
   dest->a = srce->a;
   dest->b = srce->b;
   dest->c = srce->c;
   dest->d = srce->d;
   dest->e = srce->e;
   dest->f = srce->f;
   dest->g = srce->g;
   dest->h = srce->h;
   dest->i = srce->i;
   dest->j = srce->j;
   dest->position.x = srce->position.x;
   dest->position.y = srce->position.y;
   dest->position.z = srce->position.z;
   return(dest);
}

static void primitive(code)
int code;
{
  SHAPE *shape;
  double p1, p2, p3;

  shape = alloc_shape();
  check_token(SYMLPAR);

  switch (code) {

  case SYMCUBE:
  case SYMSPHERE:
    finput(&p1, SYMRPAR, FALSE);
    p3 = p2 = p1;
    break;

  case SYMRECTANGLE:
  case SYMPLANE:
  case SYMELLIPSOID:
  case SYMCYLINDER:
  case SYMHYPERBOLOID:
  case SYMPARABOLOID:
  case SYMCONE:
    finput(&p1, SYMCOMMA, FALSE);
    finput(&p2, SYMCOMMA, FALSE);
    finput(&p3, SYMRPAR,  FALSE);
    break;

  case SYMSUR2ND:
    shape->func = _2nd_ORDER;
    finput(&(shape->a), SYMCOMMA, FALSE);
    finput(&(shape->b), SYMCOMMA, FALSE);
    finput(&(shape->c), SYMCOMMA, FALSE);
    finput(&(shape->d), SYMCOMMA, FALSE);
    finput(&(shape->e), SYMCOMMA, FALSE);
    finput(&(shape->f), SYMCOMMA, FALSE);
    finput(&(shape->g), SYMCOMMA, FALSE);
    finput(&(shape->h), SYMCOMMA, FALSE);
    finput(&(shape->i), SYMCOMMA, FALSE);
    finput(&(shape->j), SYMRPAR,  FALSE);
    break;

  default:
    error(1, IN_PRIMNAME, word);
  }

  switch (code) {

  case SYMPLANE:
    {
      vector v;

      shape->func = PLANE;
      v.x = p1;
      v.y = p2;
      v.z = p3;
      normalize(&v);
      shape->g = v.x;
      shape->h = v.y;
      shape->i = v.z;
    }
    break;

  case SYMCUBE:
  case SYMRECTANGLE:
  case SYMSPHERE:
  case SYMELLIPSOID:
  case SYMCYLINDER:
  case SYMHYPERBOLOID:
  case SYMPARABOLOID:
  case SYMCONE:
    shape->func = _2nd_ORDER;
    shape->a = sign(p1) / (p1 * p1);
    shape->b = sign(p2) / (p2 * p2);
    shape->c = sign(p3) / (p3 * p3);

    switch(code) {

    case SYMCUBE:
    case SYMRECTANGLE:
    case SYMSPHERE:
    case SYMELLIPSOID:
    case SYMCYLINDER:
    case SYMHYPERBOLOID:
      shape->j = -1.0;
      break;

    case SYMPARABOLOID:
      shape->j = 1.0;
      break;
    }
    break;
  }

  switch (code) {
  case SYMCUBE:
  case SYMRECTANGLE:
    {
      SHAPE *shsub1, *shsub2;

      shsub2 = copy_shape(alloc_shape(), shsub1 = copy_shape(alloc_shape(), shape));
      shape->b = shape->c = shsub1->a = shsub1->c = shsub2->a = shsub2->b = 0.0;
      ent_op(SYMLPAR);
      ent_operand(NEWSH_TYPE, (char *)shape);
      ent_op(OP_INTER);
      ent_operand(NEWSH_TYPE, (char *)shsub1);
      ent_op(OP_INTER);
      ent_operand(NEWSH_TYPE, (char *)shsub2);
      ent_op(SYMRPAR);
    }
    break;

  case SYMCYLINDER:
    {
      SHAPE *shsub;

      shsub = copy_shape(alloc_shape(), shape);
      shape->c = shsub->a = shsub->b = 0.0;
      ent_op(SYMLPAR);
      ent_operand(NEWSH_TYPE, (char *)shape);
      ent_op(OP_INTER);
      ent_operand(NEWSH_TYPE, (char *)shsub);
      ent_op(SYMRPAR);
    }
    break;

  default:
    ent_operand(NEWSH_TYPE, (char *)shape);
    break;
  }
}

static void enter_shape()
{
  init_expr();
  err_count = 0;

  for (;;) {
    if (toktype == OP_TOK) {
      int er;
      if (er = ent_op(opcode)) fprintf(stderr, "Parse error %d\n", er);
      if (opcode == SYMEXPTERM) break;
    }
    else if (toktype == SYMTOK) {
      if (sym == NULL) {
	union _symval temp;
	temp.objptr = (OBJECT *)NULL;
	if (ent_operand(LVALUE_TYPE, (char *)def_symbol(word, CSGSYM, &temp)))
	  error(1, BAD_SYNTAX, NULL);
      }
      else {
	switch (sym->symtype) {

	case RESSYM:
	  if (sym->symval.code >= PRCD_TOP && sym->symval.code <= PRCD_LAST)
	    primitive(sym->symval.code);
	  else if (sym->symval.code >= AFFIN_TOP && sym->symval.code <= AFFIN_LAST)
	    affin(sym->symval.code);
	  else
	    error(1, BAD_SYNTAX, NULL);
	  break;

	case CSGSYM:
	  if (ent_CSGshape(sym)) {
	    fprintf(stderr, "Illegal operand %s\n", word);
	    ++err_count;
	  }
	  break;

	case ATTSYM:
	  ent_operand(ATTR_TYPE, (char *)(sym->symval.atptr));
	  ent_op(ADDATR);
	  break;

	default:
	  fprintf(stderr, "Type mismatch of %s\n", word);
	  ++err_count;
	  break;
	}
      }
    }
    else error(1, BAD_SYNTAX, NULL);
    lexi();
  }

  switch (genCSG()) {
  case 0:
    break;

  case 1:
  case 2:
    fprintf(stderr, "Too complex\n");
    break;

  case 3:
    error(1, INVALOP, NULL);
    break;

  default:
    fprintf(stderr, "? error\n");
    break;
  }

  free_opti();
}

static void enter_sky()
{
  if (sky_def) error(1, MULT_SKY, NULL);
  sky_def = TRUE;
}

static void enter_light(symbol)
char *symbol;
{
  char token;
  double red, grn, blu;


  if (lino >= MAXNL) error(1, TOO_LIGHT, NULL);

  putc(INTLIGHT, wfp);

  switch (token = fget_token()) {

  case SYMAMBIENCE:
    putc((unsigned char)AMBI_LIT, wfp);
    break;

  case SYMPARALLEL:
    check_token(SYMLPAR);
    {
      double l1, l2;
      double p[3];
      matrix m;

      putc(PARA_LIT, wfp);
      finput(&l1, SYMCOMMA, FALSE);
      finput(&l2, SYMRPAR, FALSE);
      unit(m);
      rotx(m, DEGtoRAD(l1));
      roty(m, DEGtoRAD(l2));
      p[0] = p[1] = 0.0;
      p[2] = -1.0;
      cmat1(p, m, p);
      putd(p[0], wfp);
      putd(p[1], wfp);
      putd(p[2], wfp);
    }
    check_token(SYMCOMMA);
    break;

  case SYMPOINT:
  case SYMFPOINT:
    putc(token == SYMPOINT ? POI_LIT : FPOI_LIT, wfp);
    check_token(SYMLPAR);
    {
      double lx, ly, lz;

      finput(&lx, SYMCOMMA, FALSE);
      finput(&ly, SYMCOMMA, FALSE);
      finput(&lz, SYMRPAR, FALSE);
      putd(lx, wfp);
      putd(ly, wfp);
      putd(lz, wfp);
    }
    check_token(SYMCOMMA);
    break;

  case SYMSPOT:
    putc(SPOT_LIT, wfp);
    check_token(SYMLPAR);
    {
      double lx, ly, lz, lr, ls;
      double l1, l2, irr;
      double p[3];
      matrix m;

      finput(&lx, SYMCOMMA, FALSE);
      finput(&ly, SYMCOMMA, FALSE);
      finput(&lz, SYMCOMMA, FALSE);
      finput(&l1, SYMCOMMA, FALSE);
      finput(&l2, SYMCOMMA, FALSE);
      finput(&lr, SYMCOMMA, FALSE);
      if (lr <= 0.0) error(1, BAD_SYNTAX, NULL);
      finput(&ls, SYMRPAR, FALSE);

      irr = 1.0 / (lr * lr);
      p[0] = p[1] = 0.0;
      p[2] = -1.0;
      unit(m);
      rotx(m, DEGtoRAD(l1));
      roty(m, DEGtoRAD(l2));
      p[0] = p[1] = 0.0;
      p[2] = -1.0;
      cmat1(p, m, p);
      putd(p[0], wfp);
      putd(p[1], wfp);
      putd(p[2], wfp);
      putd(lx, wfp);
      putd(ly, wfp);
      putd(lz, wfp);
      putd(lr, wfp);
      putd(ls, wfp);
    }
    check_token(SYMCOMMA);
    break;

  default:
    error(1, BAD_SYNTAX, NULL);
  }

  red = grn = blu = 0.0;
  finput(&red, SYMCOMMA, TRUE);
  finput(&grn, SYMCOMMA, TRUE);
  finput(&blu, SYMEXPTERM, TRUE);
  putd(red, wfp);
  putd(grn, wfp);
  putd(blu, wfp);
  ++lino;
}

static int enter_camera(symbol)
char *symbol;
{
  double
    camx = 0.0, camy = 0.0, camz = 1000.0,
    tilt = 0.0, pan = 0.0, rotat = 0.0,
    zoom = 1.0;

  matrix r;
  register int i, j;


  if (camno >= 1) error(1, MULT_CAM, NULL);

  finput(&camx,  SYMCOMMA,   TRUE);
  finput(&camy,  SYMCOMMA,   TRUE);
  finput(&camz,  SYMCOMMA,   TRUE);
  finput(&tilt,  SYMCOMMA,   TRUE);
  finput(&pan,   SYMCOMMA,   TRUE);
  finput(&rotat, SYMCOMMA,   TRUE);
  finput(&zoom,  SYMEXPTERM, TRUE);
  unit(r);
  rotz(r, DEGtoRAD(rotat));
  rotx(r, DEGtoRAD(tilt));
  roty(r, DEGtoRAD(pan));
  zoom *= -50.0;
  r[2][0] *= zoom;
  r[2][1] *= zoom;
  r[2][2] *= zoom;

  putc(INTCAMERA, wfp);
  putd(camx, wfp);
  putd(camy, wfp);
  putd(camz, wfp);
  putd(zoom, wfp);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      putd(r[i][j], wfp);

  return camno++;
}

#if 0

static int enter_camera(symbol)
char *symbol;
{
  double
    camx = 0.0,
    camy = 0.0,
    camz = 0.0,
    tilt = 0.0,
    pan = 0.0,
    rotate = 0.0,
    zoom = 1.0;

  bool invert = FALSE;
  int upper = 1;

  double tarx, tary, tarz, scxx, scxy, scxz, scyx, scyy, scyz, sczx, sczy, sczz,
  vector eye;
  bool posdef, dirdef, tardef, rotdef, zoomdef;
  char token;

  if (camno >= 1) error(1, MULT_CAM, NULL);

  check_token(SYMLMPAR);
  posdef = dirdef = tardef = rotdef = zoomdef = FALSE;

  for (;;) {
    if (lexi() != SYMTOK) error(1, BAD_SYNTAX, NULL);
    if (sym == NULL) error(1, BAD_SYNTAX, NULL);
    if (sym->symtype != RESSYM) error(1, BAD_SYNTAX, NULL);
    if (( token = sym->symval.code) == SYMRMPAR) break;

    switch (token) {

    case SYMPOSITION:
      check_token(OP_EQUAL);
      if (posdef) error(0, MULT_CPOS, NULL);
      finput(&camx, SYMCOMMA, FALSE);
      finput(&camy, SYMCOMMA, FALSE);
      finput(&camz, SYMEXPTERM, FALSE);
      tilt = pan = rotate = 0.0;
      posdef = TRUE;
      break;

    case SYMTARGET:
      check_token(OP_EQUAL);
      if (dirdef || tardef) error(0, MULT_CDIR, NULL);
      dirdef = FALSE;
      finput(&tarx, SYMCOMMA, FALSE);
      finput(&tary, SYMCOMMA, FALSE);
      finput(&tarz, SYMEXPTERM, FALSE);
      tilt = pan = rotate = 0.0;
      tardef = TRUE;
      break;

    case SUMUPPER:
      check_token(OP_EQUAL);
      if (updef) error(0, MULT_DEFINE, "upper");
      token = fget_token();
      if (token == SYMPLUS || token == SYMMINUS) {
	invert = (token == SYMMINUS);
	token = fget_token();
      }
      switch (token) {
      case SYMX:  up = 0; break;
      case SYMY:  up = 1; break;
      case SYMZ:  up = 2; break;
      default:    error(1, BAD_SYNTAX, NULL);
      }
      updef = TRUE;
      break;

    case SYMROTATE:
      token = fget_token();




    case SYMTILT:
      finput(&w, SYMEXPTERM, FALSE);
      tilt += w;
      break;

    case SYMPAN:
      finput(&w, SYMEXPTERM, FALSE);
      pan += w;
      break;

    case SYMROTATE:
      finput(&w, SYMEXPTERM, FALSE);
      rotate += w;
      break;

    case SYMMOVE:
      finput(&w, SYMCOMMA, FALSE);
      camx += w;
      finput(&w, SYMCOMMA, FALSE);
      camy += w;
      finput(&w, SYMEXPTERM, FALSE);
      camz += w;
      break;

    case SYMDISTANCE:
      if (distdef) error(0, MULT_DEFINE, "distance");
      finput(&distance, SYMEXPTERM, FALSE);
      distdef = TRUE;
      break;

    case SYMZOOM:
      check_token(OP_EQUAL);
      if (zoomdef) error(0, MULT_DEFINE, "zoom");
      finput(&zoom, SYMEXPTERM, FALSE);
      zoomdef = TRUE;
      break;

    default:
      error(1, BAD_SYNTAX, NULL);
    }
  }

  eye.x = tarx - camx;
  eye.y = tary - camy;
  eye.z = tarz - camz;
  if (normalize(&eye)) error(1, CANNOT_DIRECT, NULL);

  switch (up) {

  case 0:     /* +x */
    if ((w = sqrt(mt[2][1] * mt[2][1] + mt[2][2] * mt[2][2])) == 0.0) {     /* +z */
      mt[0][0] = mt[0][2] = mt[1][0] = mt[1][1] = 0.0;
      mt[1][2] = 1.0;
      mt[0][1] = mt[2][0];
    }
    else {     /* +x */
      mt[0][0] = 0.0;
      mt[0][1] = -mt[2][2] / w;
      mt[0][2] = mt[2][1] / w;
      mt[1][0] = w;
      mt[1][1] = -mt[2][0] * mt[0][2];
      mt[1][2] = mt[2][1] * mt[0][1];
    }
    break;

  case 1:     /* +y */
    if ((w = sqrt(mt[2][0] * mt[2][0] + mt[2][2] * mt[2][2])) == 0.0) {     /* -z */
      mt[0][1] = mt[0][2] = mt[1][0] = mt[1][1] = 0.0;
      mt[1][2] = -1.0;
      mt[0][0] = mt[2][1];
    }
    else {     /* +y */
      mt[0][0] = mt[2][2] / w;
      mt[0][1] = 0.0;
      mt[0][2] = -mt[2][0] / w;
      mt[1][0] = mt[2][1] * mt[0][2];
      mt[1][1] = w;
      mt[1][2] = -mt[2][1] * mt[0][0];
    }
    break;

  case 2:     /* +z */
    if ((w = sqrt(mt[2][0] * mt[2][0] + mt[2][1] * mt[2][1])) == 0.0) {     /* +y */
      mt[0][1] = mt[0][2] = mt[1][0] = mt[1][2] = 0.0;
      mt[1][1] = 1.0;
      mt[0][0] = mt[2][2];
    }
    else {     /* +z */
      mt[0][0] = -mt[2][1] / w;
      mt[0][1] = mt[2][0] / w;
      mt[0][2] = 0.0;
      mt[1][0] = -mt[2][2] * mt[0][1];
      mt[1][1] = mt[2][2] * mt[0][0];
      mt[1][2] = w;
    }
    break;
  }

  if (invert)
    for (i = 0; i < 2; i++)
      for (j = 0; j < 3; j++)
	mt[i][j] = -mt[i][j];

  unit(r);
  rotz(r, rotate);
  rotx(r, tilt);
  roty(r, pan);
  mul_matrix(r, &mt);
  zoom *= -50.0;
  r[2][0] *= zoom;
  r[2][1] *= zoom;
  r[2][2] *= zoom;

  putc(INTCAMERA, wfp);
  putd(camx, wfp);
  putd(camy, wfp);
  putd(camz, wfp);
  putd(zoom, wfp);

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      putd(r[i][j], wfp);

  return camno++;
}

#endif

static void put_cluster()
{
  ORLIST *term;
  double putx = 0.0, puty = 0.0, putz = 0.0;

  lexi();
  if (toktype != SYMTOK || sym == NULL || sym->symtype != CSGSYM)
    error(1, BAD_SYNTAX, NULL);

  term = sym->symval.objptr->or_root;

  for (;;) {
    lexi();
    if (toktype == OP_TOK) {
      if (opcode == SYMEXPTERM) break;
      else error(1, BAD_SYNTAX, NULL);
    }
    if (toktype == SYMTOK) {
      if (sym == NULL || sym->symtype != RESSYM) error(1, BAD_SYNTAX, NULL);
      switch(sym->symval.code) {
      case SYMAT:
	check_token(SYMLPAR);
	finput(&putx, SYMCOMMA, TRUE);
	finput(&puty, SYMCOMMA, TRUE);
	finput(&putz, SYMRPAR, TRUE);
	break;

      default:
	error(1, BAD_SYNTAX, NULL);
      }
    }
    else error(1, BAD_SYNTAX, NULL);
  }

  if (term) {
    putc(INTOBJECT, wfp);
    putd(putx, wfp);
    putd(puty, wfp);
    putd(putz, wfp);

    for ( ; term; term = term->next_term) {
      ANDFACT *factor;

      putc(0, wfp);
      for (factor = term->and_root; factor; factor = factor->next_factor) {
	SHAPE *sh;
	ATTR *at;

	putc(0, wfp);
	puti((int)(sh = factor->shape)->shno, wfp);
	if (sh->shno == 0) {
	  putc(sh->func, wfp);
	  putd(sh->a, wfp);
	  putd(sh->b, wfp);
	  putd(sh->c, wfp);
	  putd(sh->d, wfp);
	  putd(sh->e, wfp);
	  putd(sh->f, wfp);
	  putd(sh->g, wfp);
	  putd(sh->h, wfp);
	  putd(sh->i, wfp);
	  putd(sh->j, wfp);
	  putd(sh->position.x, wfp);
	  putd(sh->position.y, wfp);
	  putd(sh->position.z, wfp);
	  sh->shno = shno++;
	}

	puti((at = factor->attribute)->atno, wfp);
	if (at->atno == 0) {
	  putc(at->tranflag ? 1 : 0, wfp);
	  putd(at->diff_red, wfp);
	  putd(at->diff_grn, wfp);
	  putd(at->diff_blu, wfp);
	  putd(at->spec_red, wfp);
	  putd(at->spec_grn, wfp);
	  putd(at->spec_blu, wfp);
	  putd(at->tran_red, wfp);
	  putd(at->tran_grn, wfp);
	  putd(at->tran_blu, wfp);
	  putd(at->lumi_red, wfp);
	  putd(at->lumi_grn, wfp);
	  putd(at->lumi_blu, wfp);
	  putd(at->hilit, wfp);
	  putd(at->hilit_diff, wfp);
	  putd(at->refract, wfp);
	  putc(at->texmap, wfp);
	  at->atno = atno++;
	}
	putc(factor->positive ? 1 : 0, wfp);
      }
      putc(0xff, wfp);
    }
    putc(0xff, wfp);
  }
}

static void compile(source)
char *source;
{
  FILE *fp;
  unsigned long line_count;
  FILE *m_rfp;
  char *m_inpfile;
  unsigned long *m_inpline;
  char m_ch;

  if ((fp = fopen(source, "r")) == NULL) error(1, CANT_INPUT, source);
  line_count = 1;
  m_rfp = rfp;
  rfp = fp;
  m_inpfile = inpfile;
  inpfile = source;
  m_inpline = inpline;
  inpline = &line_count;
  m_ch = ch;
  getch();

  while (lexi()) {
    if (toktype == OP_TOK && opcode == SYMEXPTERM) continue;
    if (toktype != SYMTOK) error(1, BAD_SYNTAX, NULL);
    if (sym == NULL || sym->symtype != RESSYM) enter_shape();
    else {
      switch(sym->symval.code) {

      case SYMATTRIBUTE:
	fget_attribute();
	break;

      case SYMCAMERA:
	enter_camera(NULL);
	break;

      case SYMLIGHT:
	enter_light(NULL);
	break;

      case SYMSKY:
	enter_sky ();
	break;

      case SYMPUT:
	put_cluster ();
	break;

      case SYMINCLUDE:
	{
	  char inclname[(FNAM_LEN + 1)];
	  if (lexi() != STRTOK) error(1, BAD_SYNTAX, NULL);
	  if (strlen(string) > FNAM_LEN) error(0, FNAM_TOOLONG, string);
	  compile(strncpy(inclname, string, FNAM_LEN));
	}
	break;

      default:
	error(1, BAD_SYNTAX, NULL);
      }
    }
  }
  fclose(fp);
  rfp = m_rfp;
  inpfile = m_inpfile;
  inpline = m_inpline;
  ch = m_ch;
}

static void init_object()
{
  dattrib = alloc_attr();
  dattrib->atno = 0;
  dattrib->tranflag = FALSE;	/* transparent flag */
  dattrib->diff_red   = 0.7;	/* diffusive reflection coefficient */
  dattrib->diff_grn   = 0.7;
  dattrib->diff_blu   = 0.7;
  dattrib->spec_red   = 0.001;	/* specular reflection coefficient */
  dattrib->spec_grn   = 0.001;
  dattrib->spec_blu   = 0.001;
  dattrib->tran_red   = 1.0;	/* penetration coefficient */
  dattrib->tran_grn   = 1.0;
  dattrib->tran_blu   = 1.0;
  dattrib->lumi_red   = 0.0;	/* luminous intensity */
  dattrib->lumi_grn   = 0.0;
  dattrib->lumi_blu   = 0.0;
  dattrib->hilit      = 500.0;	/* highlight coefficient */
  dattrib->hilit_diff = 0.2;	/* highlight diffusibility */
  dattrib->refract    = 1.0;	/* refractive index */
  dattrib->texmap     = 0;

  atno = shno = 1;
}

static void iscom(inpname, outname)
char *inpname, *outname;
{
  static char compiler[] = "- ISCOM -";

  camno = lino = 0;
  sky_def = FALSE;
  init_sym();
  alloc_ctree();
  init_object();
  err_count = 0;

  inpfile = compiler;
  inpline = NULL;
  outfile = outname;
  if ((wfp = fopen(outname, "w")) == NULL)
    error(1, CANT_CREATE, outname);
  compile(inpname);
  fclose(wfp);
}

int main(argc, argv)
int argc;
char *argv[];
{
  char *program = argv[0];
  bool outdes = FALSE;
  char inpname[( FNAM_LEN + 1)], outname[( FNAM_LEN + 1)];


  if (argc <= 1) {
  usage:
    fprintf(stderr, "Usage: %s [ -o output ] files\n", program);
    exit(1);
  }

  while (--argc > 0 && **++argv == '-') {
    char *flagptr = *argv;

    while (*++flagptr) {
      switch(*flagptr) {
      case 'o':
	if (outdes) {
	  fprintf(stderr, "Multiple output designation flag\n");
	  exit(1);
	}
	if (--argc == 0) {
	  fprintf(stderr, "Missing description of output\n");
	  exit(1);
	}
	strcpy(outname, *++argv);
	outdes = TRUE;
	break;

      default:
	fprintf(stderr, "Invalid flag %c\n", *flagptr);
	exit(1);
      }
    }
  }

  if (argc == 0) goto usage;

  strcpy(inpname, *argv);
  if (argc != 1) {
    fprintf(stderr, "%s may be not read\n", *++argv);
    exit(1);
  }
  if (!outdes) strcat(strcpy(outname, inpname), ".pis");

  iscom(inpname, outname);
}

#if 0

int atouie(s, termc, nextptr)
char *s, termc, *(*nextptr);
{
   register int val = 0;

   if (*s == termc)
      val = -1;
   else
   {
      while (isdigit(*s))
      {
         if (( val = 10 * val + *s++ - '0') < 0)
         {
	    val = -1;
	    break;
	 }
      }
      if (*s != termc)
	 val = -1;
   }
   *nextptr = ++s;
   return(val);
}

void printno(n, s)
unsigned int n;
char *s;
{
  if (n == 0) printf("        No");
  else printf("%10u", n);
  if (s) {
    printf(" %s", s);
    if (n > 1) putchar('s');
  }
}

#endif

/* end */
