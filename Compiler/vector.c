/* vector.c */

#include <math.h>
#include "config.h"
#include "vector.h"

double vlen(v)
vector *v;
{
  return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

bool normalize(v)
vector *v;
{
  double l;

  if (( l = vlen(v)) != 0.0) {
    v->x /= l;
    v->y /= l;
    v->z /= l;
    return FALSE;
  }
  else return TRUE;
}

/* end */
