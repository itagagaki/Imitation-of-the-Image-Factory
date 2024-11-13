#include <stdio.h>
#include <math.h>
#include "config.h"
#include "vector.h"
#include "object.h"
#include "matrix.h"
#include "error.h"

static matrix temp;

matrixp alloc_matrix()
{
  matrixp m;

  if ((m = (matrixp)malloc(sizeof(matrix))) == NULL)
    error(1, TOO_SHAPE, NULL);

  return m;
}

matrixp unit(m)
matrixp m;
{
  register int j, i;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      m[i][j] = j == i ? 1.0 : 0.0;

  return m;
}

matrixp mul_matrix(m1, m2)
matrixp m1, m2;
{
  register int j, i;
  double a, b, c, d;

  for (i = 0; i < 4; i++) {
    a = m1[i][0];
    b = m1[i][1];
    c = m1[i][2];
    d = m1[i][3];
    for (j = 0; j < 4; j++)
      m1[i][j] = a * m2[0][j] + b * m2[1][j] + c * m2[2][j] + d * m2[3][j];
  }

  return m1;
}

matrixp rotxp(m, sinr, cosr)
matrixp m;
double sinr, cosr;
{
  unit(temp);
  temp[1][1] = temp[2][2] = cosr;
  temp[2][1] = -(temp[1][2] = sinr);
  return mul_matrix(m, temp);
}

matrixp rotyp(m, sinr, cosr)
matrixp m;
double sinr, cosr;
{
  unit(temp);
  temp[0][0] = temp[2][2] = cosr;
  temp[0][2] = -(temp[2][0] = sinr);
  return mul_matrix(m, temp);
}

matrixp rotzp(m, sinr, cosr)
matrixp m;
double sinr, cosr;
{
  unit(temp);
  temp[0][0] = temp[1][1] = cosr;
  temp[1][0] = -(temp[0][1] = sinr);
  return mul_matrix(m, temp);
}

matrixp rotx(m, r)
matrixp m;
double r;
{
  return rotxp(m, sin(r), cos(r));
}

matrixp roty(m, r)
matrixp m;
double r;
{
  return rotyp(m, sin(r), cos(r));
}

matrixp rotz(m, r)
matrixp m;
double r;
{
  return rotzp(m, sin(r), cos(r));
}

matrixp move(m, mv)
matrixp m;
vector *mv;
{
   unit(temp);
   temp[3][0] = mv->x;
   temp[3][1] = mv->y;
   temp[3][2] = mv->z;
   return mul_matrix(m, temp);
}

matrixp scale(m, sv)
matrixp m;
vector *sv;
{
  unit(temp);
  temp[0][0] = sv->x;
  temp[1][1] = sv->y;
  temp[2][2] = sv->z;
  return mul_matrix(m, temp);
}

double *cmat1(sp, m, dp)
double sp[3], dp[3];
matrixp m;
{
  double a, b, c;

  a = sp[0];
  b = sp[1];
  c = sp[2];
  dp[0] = a * m[0][0] + b * m[1][0] + c * m[2][0];
  dp[1] = a * m[0][1] + b * m[1][1] + c * m[2][1];
  dp[2] = a * m[0][2] + b * m[1][2] + c * m[2][2];

  return dp;
}

#if 0

double *cmat2(sp, m, dp)
double sp[6], dp[6];
matrixp m;
{
  double a, b, c, d, e, f;

  a = sp[0];
  b = sp[1];
  c = sp[2];
  d = sp[3];
  e = sp[4];
  f = sp[5];

  dp[0] = a * m[0][0] * m[0][0]
        + b * m[1][0] * m[1][0]
	+ c * m[2][0] * m[2][0]
	+ d * m[1][0] * m[2][0]
	+ e * m[0][0] * m[2][0]
	+ f * m[0][0] * m[1][0];

  dp[1] = a * m[0][1] * m[0][1]
        + b * m[1][1] * m[1][1]
	+ c * m[2][1] * m[2][1]
	+ d * m[1][1] * m[2][1]
	+ e * m[0][1] * m[2][1]
	+ f * m[0][1] * m[1][1];

  dp[2] = a * m[0][2] * m[0][2]
        + b * m[1][2] * m[1][2]
	+ c * m[2][2] * m[2][2]
	+ d * m[1][2] * m[2][2]
	+ e * m[0][2] * m[2][2]
	+ f * m[0][2] * m[1][2];

  dp[3] = 2.0 * (a * m[0][1] * m[0][2]
	       + b * m[1][1] * m[1][2]
	       + c * m[2][1] * m[2][2])
        + d * (m[1][1] * m[2][2] + m[2][1] * m[1][2])
	+ e * (m[0][1] * m[2][2] + m[2][1] * m[0][2])
	+ f * (m[0][1] * m[1][2] + m[1][1] * m[0][2]);

  dp[4] = 2.0 * (a * m[0][0] * m[0][2]
	       + b * m[1][0] * m[1][2]
	       + c * m[2][0] * m[2][2])
        + d * (m[1][0] * m[2][2] + m[2][0] * m[1][2])
	+ e * (m[0][0] * m[2][2] + m[2][0] * m[0][2])
	+ f * (m[0][0] * m[1][2] + m[1][0] * m[0][2]);

  dp[5] = 2.0 * (a * m[0][0] * m[0][1]
	       + b * m[1][0] * m[1][1]
	       + c * m[2][0] * m[2][1])
        + d * (m[1][0] * m[2][1] + m[2][0] * m[1][1])
	+ e * (m[0][0] * m[2][1] + m[2][0] * m[0][1])
	+ f * (m[0][0] * m[1][1] + m[1][0] * m[0][1]);

  return dp;
}

#endif

vector *apply(v, m)
vector *v;
matrixp m;
{
  double x, y, z, w;

  x = v->x * m[0][0] + v->y * m[1][0] + v->z * m[2][0] + m[3][0];
  y = v->x * m[0][1] + v->y * m[1][1] + v->z * m[2][1] + m[3][1];
  z = v->x * m[0][2] + v->y * m[1][2] + v->z * m[2][2] + m[3][2];
  w = v->x * m[0][3] + v->y * m[1][3] + v->z * m[2][3] + m[3][3];
  v->x = x / w;
  v->y = y / w;
  v->z = z / w;

  return v;
}

/* end */
