#include <stdio.h>
#include <math.h>
#include "config.h"
#include "form.h"
#include "vector.h"
#include "object.h"

unsigned char cred[SXMAX], cgrn[SXMAX], cblu[SXMAX], cmask[SXMAX];

extern int super_level, super_n;
extern double camera_x, camera_y, camera_z;
extern double cos_camtilt, sin_camtilt, cos_campan, sin_campan, cos_camrotate, sin_camrotate;
extern double focal_distance;
extern double vpx, vpy, vpz;
extern int reso_h;
extern double screen_left_edge, screen_upper_edge;
extern double pixel_size_x, pixel_size_y, pixel_size_xs, pixel_size_ys;

extern double total_red, total_grn, total_blu;

extern bool trace();


/*
	compute the shade value for each pixel of a laster
*/
void compute(laster)
int laster;
{
  register int j, i;
  double yo1sv3[MAXSUPER];
  double yo1cv3[MAXSUPER];
  double met1[MAXSUPER];
  double yss12[MAXSUPER];
  double ysc12[MAXSUPER];
  double xo1;
  double metric;
  double vx, vy, vz;
  double temp;
  unsigned char *pred, *pgrn, *pblu;

  xo1 = screen_left_edge;
  pred = cred;
  pgrn = cgrn;
  pblu = cblu;

  lox = 

  xo1 = screen_left_edge + pixel_size_xs / 2.0;
  i = reso_h;
  while (i--) {
    total_red = total_grn = total_blu = 0.0;
    for (j = 0; j < super_level; j++) {
      register int k;

      for (k = 0; k < super_level; k++) {
	double l;

	vx = sox + (svx - sox) / reso_v * laster + (svx - sox) / reso_v / (super_level + 1) * (k + 1) - ex;
	vx = sox + ((svx - sox) / reso_v) * (laster + (k + 1) / (super_level + 1)) - ex;
	vy = py - ey;
	vz = pz - ez;
	l = sqrt(vx * vx + vy * vy + vz * vz);
	trace(camera_x, camera_y, camera_z,
	      vx / l, vy / l, vz / l,
	      NULL,
	      1,
	      1.0, 1.0, 1.0,
	      1.0,
	      0.0, 0.0, 0.0);
      }
      xo1 += pixel_size_xs;
    }
    total_red /= super_n;
    total_grn /= super_n;
    total_blu /= super_n;
    *pred++ = total_red <= 0.0 ? 0 : (total_red >= 1.0 ? 255 : total_red * 256.0);
    *pgrn++ = total_grn <= 0.0 ? 0 : (total_grn >= 1.0 ? 255 : total_grn * 256.0);
    *pblu++ = total_blu <= 0.0 ? 0 : (total_blu >= 1.0 ? 255 : total_blu * 256.0);
  }
}

#if 0

double marume(value, base)
double value, base;
{
  double in;

  if (modf((value / base), &in) < 0.5) return in * base;
  else return (in + (in > 0.0 ? 1.0 : -1.0)) * base;
}

#endif

#if 0
pp(ptr)
unsigned char *ptr;
{
  int i, j;
  static unsigned int mb[8] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};

  ptr += 7;

  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) putchar((*ptr & mb[j]) == 0 ? '0' : '1');
    --ptr;
  }

  putchar('\n');
}

#endif
