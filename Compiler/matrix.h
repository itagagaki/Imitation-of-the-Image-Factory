typedef double matrix[4][4];
typedef double (*matrixp)[4];

extern matrixp
  alloc_matrix(),
  unit(),
  mul_matrix(),
  rotxp(),
  rotyp(),
  rotzp(),
  rotx(),
  roty(),
  rotz(),
  move(),
  scale();
extern double *cmat1(), *cmat2();
extern vector *apply();
