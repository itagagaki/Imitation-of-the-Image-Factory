#include <stdio.h>
#include <vmath.h>
#include "config.h"
#include "form.h"
#include "math2.h"
#include "vector.h"
#include "object.h"
#include "light.h"
#include "toleran.h"

extern FILE *wfp;
unsigned char scan_line_buffer[SXMAX*4];

extern int anti_depth;
extern int anti_thres;
extern int obj_complex_flag;
extern int cam_in_obj_flag;
extern int super_level;
extern double camera_x, camera_y, camera_z, focal_dist,
	      ca11, ca12, ca13, ca21, ca22, ca23,
	      ca31z, ca32z, ca33z;

extern int reso_h, reso_v;
extern double screen_left_edge, screen_upper_edge;
extern double pixel_size_x, pixel_size_y;

double cam_n, cam_ar, cam_ag, cam_ab;
int sample_red, sample_grn, sample_blu;
double sub_red, sub_grn, sub_blu;
double total_red, total_grn, total_blu;


extern double min_expseed;
extern int shadow_sw, hilit_sw;
extern int trace_count;
extern double min_energy;
extern unsigned int num_light;

extern char light_kind[MAXNL];
extern double light_x[MAXNL],light_y[MAXNL],light_z[MAXNL],
	      light_vx[MAXNL],light_vy[MAXNL],light_vz[MAXNL],
	      light_r[MAXNL], light_s[MAXNL],
	      light_red[MAXNL],light_grn[MAXNL],light_blu[MAXNL];

extern CLUSTER *object_root;

double x0_objloc, y0_objloc, z0_objloc;
double pre_x1_objloc, pre_y1_objloc, pre_z1_objloc;
double diff_red, diff_grn, diff_blu;

double ray_red, ray_grn, ray_blu;


#define sfabs(d) ((d)<0.0?-(d):(d))

#define sq(d) ((d)*(d))

#define vabs(x,y,z) sqrt(sq(x)+sq(y)+sq(z))

#define normalize(w,x,y,z) if ((w = vabs(x,y,z)) != 0.0) ((x/=w),(y/=w),(z/=w))

/*
	specificate the color at intersection
*/
void texture(attr, x, y, z)
ATTR *attr;
double x, y, z;
{
  diff_red = attr->diff_red;
  diff_grn = attr->diff_grn;
  diff_blu = attr->diff_blu;

  switch (attr->texmap)
   {
   case 1:
     if ( ((x - floor(x * 0.05) * 20.0) < 10.0)
	^ ((y - floor(y * 0.05) * 20.0) < 10.0)
	^ ((z - floor(z * 0.05) * 20.0) < 10.0) )
       diff_grn = 0.0;
     break;

   case 2:
     if ( ((x - floor(x / 60.0) * 60.0) < 30.0)
	^ ((y - floor(y / 60.0) * 60.0) < 30.0)
	^ ((z - floor(z / 60.0) * 60.0) < 30.0) ) {
       diff_red = 0.4;
       diff_grn = 0.55;
     }
     break;

   case 3:
     if ( ((x - floor(x * 0.05) * 20.0) < 10.0)
	^(( y - floor(y * 0.05) * 20.0) < 10.0)
	^(( z - floor(z * 0.05) * 20.0) < 10.0) )
       diff_grn = 0.0;
     break;

   case 4:
     if ( ((x - floor(x / 60.0) * 60.0) < 30.0)
	^ (( y - floor(y / 60.0) * 60.0) < 30.0)
	^ (( z - floor(z / 60.0) * 60.0) < 30.0) ) {
       diff_red *= 0.4;
       diff_grn *= 0.55;
     }
     break;

   case 0:
   default:
     break;
   }
}

/*
	is the point out side of the object?
*/
bool isout(primitive)
ANDFACT *primitive;
{
  double q;
  SHAPE *sh;
  double x, y, z;

  sh = primitive->shape;
  x = pre_x1_objloc - sh->position.x;
  y = pre_y1_objloc - sh->position.y;
  z = pre_z1_objloc - sh->position.z;
  q = sh->a * sq(x) + sh->b * sq(y) + sh->c * sq(z)
    + sh->d * y * z + sh->e * x * z + sh->f * x * y
    + sh->g * x + sh->h * y + sh->i * z
    + sh->j;
  return primitive->positive ? q > 0.0 : q <= 0.0;
}

/*
	calculate inter section
*/
int intsct(primitive, vx, vy, vz, dist)
ANDFACT *primitive;
double vx, vy, vz;
double dist[];
{
/*
ray equation:
     x = vx*t + x0
     y = vy*t + y0
     z = vz*t + z0
*/

  SHAPE *cur_shape;
  double qa, qb, qc, qd;
  double x, y, z;

  cur_shape = primitive->shape;
  x = x0_objloc - cur_shape->position.x;
  y = y0_objloc - cur_shape->position.y;
  z = z0_objloc - cur_shape->position.z;

  switch (cur_shape->func)
    {
    case PLANE:
/*
surface equation:
     a*x + b*y + c*z = 0
*/
      if ((qd = vx * cur_shape->g + vy * cur_shape->h + vz * cur_shape->i) == 0.0)
	return 0;
      else {
	dist[0] = -(x * cur_shape->g + y * cur_shape->h + z * cur_shape->i) / qd;
	return 1;
      }

    case _2nd_ORDER:
/*
surface equation:
   a*x^2 + b*y^2 + c*z^2 + d*y*z + e*x*z + f*x*y + g*x + h*y + i*z + j = 0

simultaneous ray equation and this:
   qa*t^2 + 2*qb*t + qc = 0

   qa = a*vx^2 + b*vy^2 + c*vz^2 + d*vy*vz + e*vx*vz + f*vx*vy

   qb = a*vx*x0 + b*vy*y0 + c*vz*z0 +(d*(vy*z0 + vz*y0) + e*(vxz0 + vzx0) + f*(vxy0 + vyx0) + g*vx + h*vy + i*vz)/2

   qc = a*x0^2 + b*y0^2 + c*z0^2 + d*y0*z0 + e*x0*z0 + f*x0*y0 + j

discriminant:	qd = qb^2 - qa * qc

solution:
		t1 = -qb' - sqrt(qd) / qa

		t2 = -qb' + sqrt(qd) / qa

		t1*t2 = qc / qa
*/
      qa = cur_shape->a * sq(vx)
	 + cur_shape->b * sq(vy)
	 + cur_shape->c * sq(vz)
	 + cur_shape->d * vy * vz
	 + cur_shape->e * vx * vz
	 + cur_shape->f * vx * vy;

      if (qa == 0.0) return 0;

      qb = cur_shape->a * vx * x
	 + cur_shape->b * vy * y
	 + cur_shape->c * vz * z
	 + ( cur_shape->d * (z * vy + y * vz)
	   + cur_shape->e * (x * vz + z * vx)
	   + cur_shape->f * (x * vy + y * vx)
	   + cur_shape->g * vx
	   + cur_shape->h * vy
	   + cur_shape->i * vz
	   ) * 0.5;

      qc = cur_shape->a * sq(x)
	 + cur_shape->b * sq(y)
	 + cur_shape->c * sq(z)
         + cur_shape->d * y * z
	 + cur_shape->e * x * z
	 + cur_shape->f * x * y
	 + cur_shape->j;

      if ((qd = sq(qb) - qa * qc) >= 0.0) {
	if ( (dist[0] =
	        (qd = (qb > 0.0 ? -sqrt(qd) : sqrt(qd)) - qb) / qa) == 0.0 )
	  return 1;

	dist[1] = qc / qd;
	if (dist[0] > dist[1]) {
	  double tmp = dist[0];
	  dist[0] = dist[1];
	  dist[1] = tmp;
	}
	return 2;
      }
      else return 0;

    default:
      return 0;

    }
}

/*
	search for object hits the ray
*/
bool search(x0, y0, z0, vx, vy, vz, finite, d_lim,
	    last_factor, last_enter,
	    distance, x1, y1, z1, x1f, y1f, z1f,
	    pofx, pofy, pofz,
	    object, factor, enter, right)
double x0, y0, z0;
double vx, vy, vz;
bool finite;
double d_lim;
ANDFACT *last_factor;
bool last_enter;
double *distance;
double *x1, *y1, *z1;
double *x1f, *y1f, *z1f;
double *pofx, *pofy, *pofz;
CLUSTER **object;
ANDFACT **factor;
bool *enter, *right;
{
  bool found = FALSE;
  CLUSTER *objptr, *cur_object;
  ORLIST *orptr;

  ANDFACT *cur_factor;
  double tmin;
  bool cur_enter;

  for (objptr = object_root; objptr; objptr = objptr->next_cluster) {
    x0_objloc = x0 - objptr->putx;
    y0_objloc = y0 - objptr->puty;
    z0_objloc = z0 - objptr->putz;

#if 0
    if (range_prim != 999) {
      if (intsct()) continue;
      if (distance >= tmin) continue;
    }
#endif

    for (orptr = objptr->or_root; orptr; orptr = orptr->next_term) {
      ANDFACT *and_root, *andptr;

      for (andptr = and_root = orptr->and_root; andptr; andptr = andptr->next_factor) {
	int num_intsct;
	double dist[2];
	register int chklc;

	if ((num_intsct = intsct(andptr, vx, vy, vz, dist)) < 1) {
#if 0
	  if (andptr->positive) break;
	  else continue;
#else
	  continue;
#endif
	}

	/* Fail if the smallest is far from decided one. */
	if (found) {
	  if (dist[0] - NEARTOL > tmin) continue;
	}

	/* Fail if the smallest is far from the limit. */
	if (finite) {
	  if (dist[0] - NEARTOL > d_lim) continue;
	}

	/* Fail if the largest is back of the starting point. */
	if (dist[num_intsct - 1] < -ZEROTOL) continue;

	/* Cut the nearest if this is last primitive. */
	if ((char *)andptr == (char *)last_factor) {
	  double am;
	  int id, is;

	  id = 0;
	  am = sfabs(dist[0]);

	  for (is = 1; is < num_intsct; is++)
	    if (sfabs(dist[is]) < am) {
	      id = is;
	      am = sfabs(dist[is]);
	    }

	  is = id + 1;
	  while (is < num_intsct) dist[id++] = dist[is++];
	  --num_intsct;
	}

	/* Test each intersection. */

	for (chklc = 0; chklc < num_intsct; chklc++) {
	  double di, dii;
	  bool enter;

	  /* Fail if intersection is back of the starting point. */
	  if ((di = dist[chklc]) < -ZEROTOL) continue;

	  /* Fail if intersection is just starting point
	     and the point was already touched until now. */
	  if (di < ZEROTOL) {
	    if (last_enter) continue;
	  }

	  /* Fail if the smallest is far from limit. */
	  if (finite) {
	    if (di - NEARTOL > d_lim) continue;
	  }

	  /* Fail if intersection is from decided one. */
	  if (found) {
	    if (di - NEARTOL > tmin) continue;
	  }

	  dii = di + MARGIN;
	  pre_x1_objloc = x0_objloc + vx * dii;
	  pre_y1_objloc = y0_objloc + vy * dii;
	  pre_z1_objloc = z0_objloc + vz * dii;
	  enter = !isout(andptr);

	  if (found) {
	    if (di + NEARTOL > tmin) {
	      if (!cur_factor->attribute->tranflag) continue;
	      if (andptr->attribute->tranflag) {
		if (cur_enter) continue;
		if (!enter) continue;
	      }
	    }
	  }

	  if (!enter) {
	    if (di < ZEROTOL) continue;
	    dii = di - MARGIN;
	    pre_x1_objloc = x0_objloc + vx * dii;
	    pre_y1_objloc = y0_objloc + vy * dii;
	    pre_z1_objloc = z0_objloc + vz * dii;
	  }

	  /* test AND operation */

	  {
	    ANDFACT *andref;
	    register bool not_fail = TRUE;

	    for (andref = and_root; andref; andref = andref->next_factor) {
	      if (andref != andptr) {
		if (isout(andref)) {
		  not_fail = FALSE;
		  break;
		}
	      }
	    }

	    if (not_fail) {
	      tmin = di;
	      cur_enter = enter;
	      cur_object = objptr;
	      cur_factor = andptr;
	      found = TRUE;
	      break;
	    }
	  }
	}
      }
    }
  }

  if (found) {
    double dif;

    *distance = tmin < 0.0 ? 0.0 : tmin;
    *x1 = x0 + vx * tmin;
    *y1 = y0 + vy * tmin;
    *z1 = z0 + vz * tmin;
    dif = tmin + MARGIN;
    *x1f = x0 + vx * dif;
    *y1f = y0 + vy * dif;
    *z1f = z0 + vz * dif;
    *pofx =(*object = cur_object)->putx + (*factor = cur_factor)->shape->position.x;
    *pofy = cur_object->puty + cur_factor->shape->position.y;
    *pofz = cur_object->putz + cur_factor->shape->position.z;
    *enter = cur_enter;
    *right = cur_factor->positive ? cur_enter : !cur_enter;
  }

  return found;
}

/*
	calc complex refract
*/
void calc_refract(object, x1, y1, z1, refract, red, grn, blu)
CLUSTER *object;
double x1, y1, z1;
double *refract, *red, *grn, *blu;
{
  ORLIST *orptr;
  ANDFACT *andptr;
  ATTR *at;
  int denomi;

  *refract = *red = *grn = *blu = 0.0;
  denomi = 0;
  pre_x1_objloc = x1 - object->putx;
  pre_y1_objloc = y1 - object->puty;
  pre_z1_objloc = z1 - object->putz;
  for (orptr = object->or_root; orptr; orptr = orptr->next_term) {
    double sub_refr, sub_ar, sub_ag, sub_ab;
    int sub_deno;

    sub_refr = sub_ar = sub_ag = sub_ab = 0.0;
    sub_deno = 0;

    for (andptr = orptr->and_root;
	 andptr && !isout(andptr);
	 andptr = andptr->next_factor)
      if ((at = andptr->attribute)->tranflag) {
	sub_refr += at->refract;
	sub_ar += at->tran_red;
	sub_ag += at->tran_grn;
	sub_ab += at->tran_blu;
	++sub_deno;
      }

    if (andptr == NULL) {
      *refract += sub_refr;
      *red += sub_ar;
      *grn += sub_ag;
      *blu += sub_ab;
      denomi += sub_deno;
    }
  }

  if (denomi == 0) {
    *refract = 1.0;
    *red = *grn = *blu = 0.0;
  }
  else {
    *refract /= denomi;
    *red /= denomi;
    *grn /= denomi;
    *blu /= denomi;
  }
}

bool whereis( x1, y1, z1, vx, vy, vz, factor, out_obj, out_fact)
double x1, y1, z1;
double vx, vy, vz;
ANDFACT *factor;
CLUSTER **out_obj;
ANDFACT **out_fact;
{
  CLUSTER *object;
  double x2, y2, z2;
  ANDFACT *factor2;
  bool enter, last_enter;
  unsigned int count;
  double d_dist;
  double d_x1f, d_y1f, d_z1f, d_pofx, d_pofy, d_pofz;
  bool d_rev;

  last_enter = FALSE;
  count = 1;
  do {
    if (search(x1, y1, z1, vx, vy, vz, FALSE, 0.0,
	       factor, last_enter,
	       &d_dist, &x2, &y2, &z2, &d_x1f, &d_y1f, &d_z1f,
	       &d_pofx, &d_pofy, &d_pofz,
	       &object, &factor2, &enter, &d_rev)) {
      if (enter) ++count;
      else --count;
      x1 = x2;
      y1 = y2;
      z1 = z2;
      factor = factor2;
      last_enter = enter;
    }
    else return FALSE;
  } while (count > 0);

  *out_obj = object;
  *out_fact = factor2;
  return TRUE;
}

/*
	calculate normal vector at intersection
*/
void nvector(shape, posx, posy, posz, right, nx, ny, nz)
SHAPE *shape;
double posx, posy, posz;
bool right;
double *nx, *ny, *nz;
{
  switch (shape->func)
    {
    case PLANE:
      if (right) {
	*nx = shape->g;
	*ny = shape->h;
	*nz = shape->i;
      }
      else {
	*nx = -shape->g;
	*ny = -shape->h;
	*nz = -shape->i;
      }
      break;

    case _2nd_ORDER:
      {
	double w, vx, vy, vz;

	vx = shape->a * posx + (shape->e * posz + shape->f * posy + shape->g) * 0.5;

	vy = shape->b * posy + (shape->d * posz + shape->f * posx + shape->h) * 0.5;

	vz = shape->c * posz + (shape->d * posy + shape->e * posx + shape->i) * 0.5;
	normalize(w, vx, vy, vz);
	if (right) {
	  *nx = vx;
	  *ny = vy;
	  *nz = vz;
	}
	else {
	  *nx = -vx;
	  *ny = -vy;
	  *nz = -vz;
	}
      }
      break;
   }
}

/*
	trace the shadow
*/
bool shadow(x0, y0, z0, vx, vy, vz, finlit, ld, last_factor, last_enter,
	    n1, ar1, ag1, ab1, level, lred, lgrn, lblu)
double x0, y0, z0;
double vx, vy, vz;
bool finlit;
double ld;
ANDFACT *last_factor;
bool last_enter;
double n1;
double ar1, ag1, ab1;
int level;
double *lred, *lgrn, *lblu;
{
  double distance;
  double x1, y1, z1, x1f, y1f, z1f, pofx, pofy, pofz;
  CLUSTER *object;
  ANDFACT *factor;
  bool enter, right;
  double nx, ny, nz;
  double cosvn;
  double n2, nr;
  double ar2, ag2, ab2;
  double w;
  ATTR *at;

  if (search(x0, y0, z0, vx, vy, vz, finlit, ld,
	     last_factor, last_enter,
	     &distance,
	     &x1, &y1, &z1, &x1f, &y1f, &z1f, &pofx, &pofy, &pofz,
	     &object, &factor, &enter, &right))
    {
      if (factor->attribute->tranflag) {
	if (level == 0) {
	  CLUSTER *out_object;
	  ANDFACT *out_factor;

	  if (whereis(x0, y0, z0, vx, vy, vz, last_factor,
		      &out_object, &out_factor)) {
	    if (obj_complex_flag)
	      calc_refract(out_object,
			   x0 + vx * MARGIN,
			   y0 + vy * MARGIN,
			   z0 + vz * MARGIN,
			   &n1, &ar1, &ag1, &ab1);
	    else {
	      n1 = (at = out_factor->attribute)->refract;
	      ar1 = at->tran_red;
	      ag1 = at->tran_grn;
	      ab1 = at->tran_blu;
	    }
	  }
	  else {
	    n1 = 1.0;
	    ar1 = ag1 = ab1 = 0.0;
	  }
	}

	*lred *= (w = ar1 * -distance) >= min_expseed ? exp(w) : 0.0;
	*lgrn *= (w = ag1 * -distance) >= min_expseed ? exp(w) : 0.0;
	*lblu *= (w = ab1 * -distance) >= min_expseed ? exp(w) : 0.0;

	if (enter) {
	  if (obj_complex_flag)
	    calc_refract(object, x1f, y1f, z1f, &n2, &ar2, &ag2, &ab2);
	  else {
	    n2 = (at = factor->attribute)->refract;
	    ar2 = at->tran_red;
	    ag2 = at->tran_grn;
	    ab2 = at->tran_blu;
	  }
	}
	else {
	  CLUSTER *out_object;
	  ANDFACT *out_factor;

	  if (whereis( x1, y1, z1, vx, vy, vz, factor, &out_object, &out_factor)) {
	    if (obj_complex_flag)
	      calc_refract(out_object, x1f, y1f, z1f,
			   &n2, &ar2, &ag2, &ab2);
	    else {
	      n2 = (at = out_factor->attribute)->refract;
	      ar2 = at->tran_red;
	      ag2 = at->tran_grn;
	      ab2 = at->tran_blu;
	    }
	  }
	  else {
	    n2 = 1.0;
	    ar2 = ag2 = ab2 = 0.0;
	  }
	}

	shadow(x1, y1, z1, vx, vy, vz,
	       finlit, finlit ? ld - distance : ld,
	       factor, enter,
	       n2, ar2, ag2, ab2, level + 1, lred, lgrn, lblu);
#if 0
	nvector(factor->shape,
		x1 - pofx, y1 - pofy, z1 - pofz, !right,
		&nx, &ny, &nz );
	cosvn = vx * nx + vy * ny + vz * nz;

	if (enter) {
	  if (obj_complex_flag)
	    calc_refract(object, x1f, y1f, z1f, &n2, &ar2, &ag2, &ab2);
	  else {
	    n2 = (at = factor->attribute)->refract;
	    ar2 = at->tran_red;
	    ag2 = at->tran_grn;
	    ab2 = at->tran_blu;
	  }
	}
	else {
	  CLUSTER *out_object;
	  ANDFACT *out_factor;

	  if (whereis(x1, y1, z1, vx, vy, vz, factor,
			   &out_object, &out_factor)) {
	    if (obj_complex_flag)
	      calc_refract(out_object, x1f, y1f, z1f,
			   &n2, &ar2, &ag2, &ab2);
	    else {
	      n2 = (at = out_factor->attribute)->refract;
	      ar2 = at->tran_red;
	      ag2 = at->tran_grn;
	      ab2 = at->tran_blu;
	    }
	  }
	  else {
	    n2 = 1.0;
	    ar2 = ag2 = ab2 = 0.0;
	  }
	}

	nr = n1 / n2;

	if ((w = sq(nr) + sq(cosvn) - 1.0) >= 0.0) {
	  double ncosvn, cosmnp, ncosmnp, ss, sp, kt;

	  w = sqrt(w) - cosvn;
	  ncosvn = nr * cosvn;
	  ncosmnp = nr * (cosmnp =
			  nx * (vx + nx * w) / nr
			+ ny * (vy + ny * w) / nr
			+ nz * (vz + nz * w) / nr);
	  ss = (cosvn - ncosmnp) / (cosvn + ncosmnp);
	  sp = (ncosvn - cosmnp) / (ncosvn + cosmnp);
	  kt = 1.0 - (sq(ss) + sq(sp)) * 0.5;
	  *lred *= kt;
	  *lgrn *= kt;
	  *lblu *= kt;
	  shadow(x1, y1, z1, vx, vy, vz,
		 finlit, finlit ? ld - distance : ld,
		 factor, enter,
		 n2, ar2, ag2, ab2, level + 1, lred, lgrn, lblu);
	}
	else *lred = *lgrn = *lblu = 0.0;
#endif
      }
      else *lred = *lgrn = *lblu = 0.0;

      return TRUE;
    }
  else return FALSE;
}

/*
	trace the ray
*/
int trace(x0, y0, z0, vx, vy, vz, last_factor, last_enter,
		      level, en_red, en_grn, en_blu, n1, ar1, ag1, ab1)
double x0, y0, z0;
double vx, vy, vz;			/* view vector */
ANDFACT *last_factor;
bool last_enter;
int level;				/* generation of ray */
double en_red, en_grn, en_blu;		/* energy of ray */
double n1;				/* entrance refractive index */
double ar1, ag1, ab1;
{
  double distance;
  double x1, y1, z1;
  double x1f, y1f, z1f;
  double pofx, pofy, pofz;
  double nx, ny, nz;		/* normal vector of surface */
  CLUSTER *object;
  ANDFACT *factor;
  bool enter, right;
  ATTR *attr;
  double cosvn;
  double tmp_dred, tmp_dgrn, tmp_dblu;
  double tmp_sred, tmp_sgrn, tmp_sblu;
  double w, w1, w2, w3;
  int m;

  if (!search(x0, y0, z0, vx, vy, vz, FALSE, 0.0,
	      last_factor, last_enter,
	      &distance, &x1, &y1, &z1, &x1f, &y1f, &z1f,
	      &pofx, &pofy, &pofz,
	      &object, &factor, &enter, &right))
    return 0;

  en_red *= (w = ar1 * -distance) >= min_expseed ? exp(w) : 0.0;
  en_grn *= (w = ag1 * -distance) >= min_expseed ? exp(w) : 0.0;
  en_blu *= (w = ab1 * -distance) >= min_expseed ? exp(w) : 0.0;

  nvector(factor->shape, x1 - pofx, y1 - pofy, z1 - pofz, right, &nx, &ny, &nz);
  cosvn = -vx * nx - vy * ny - vz * nz;

  attr = factor->attribute;

  if (enter) {
    texture(attr, x1f - pofx, y1f - pofy, z1f - pofz);
    tmp_dred = tmp_dgrn = tmp_dblu = tmp_sred = tmp_sgrn = tmp_sblu = 0.0;

    for (m = 0; m < num_light; m++) {
      bool finlit;		     /* distance finity flag */
      vector lv;   /* L : vector of light */
      vector hv;   /* H : half direction vector between L and V */
      double ld;		     /* distance for light */
      double lred, lgrn, lblu;    /* light color */
      double cosln, coshn, cosshn, coshv;
      double hl;

      lred = light_red[m];
      lgrn = light_grn[m];
      lblu = light_blu[m];
      switch (light_kind[m])
	{
	case AMBI_LIT:
	  tmp_dred += lred;
	  tmp_dgrn += lgrn;
	  tmp_dblu += lblu;
	  continue;

	case PARA_LIT:
	  lv.x = light_vx[m];
	  lv.y = light_vy[m];
	  lv.z = light_vz[m];
	  finlit = FALSE;
	  break;

	case POI_LIT:
	case FPOI_LIT:
	  lv.x = x1 - light_x[m];
	  lv.y = y1 - light_y[m];
	  lv.z = z1 - light_z[m];
	  normalize(ld, lv.x, lv.y, lv.z);
	  w = ld;
	  if (light_kind[m] == POI_LIT) w *= w;
	  w += 1.0;
	  lred /= w;
	  lgrn /= w;
	  lblu /= w;
	  finlit = TRUE;
	  break;

	case SPOT_LIT:
	  {
	    double x, y, z;

	    x = x1 - light_x[m];
	    y = y1 - light_y[m];
	    z = z1 - light_z[m];
	    lv.x = light_vx[m];
	    lv.y = light_vy[m];
	    lv.z = light_vz[m];

	    ld = lv.x * x + lv.y * y + lv.z * z;
	    if (ld < 0.0) continue;   /* back of light */

	    w1 = x - lv.x * ld;
	    w2 = y - lv.y * ld;
	    w3 = z - lv.z * ld;

	    w = vabs(w1, w2, w3) / light_r[m];
	    if (w > 1.0) continue;   /* out of spot area */

	    w = pow(cos(w * PI * 0.5), light_s[m]);
	    lred *= w;
	    lgrn *= w;
	    lblu *= w;
	    finlit = TRUE;
	  }
	  break;
	}

      if ((cosln = -lv.x * nx - lv.y * ny - lv.z * nz) > 0.0) {
	if (shadow_sw)
	  shadow(x1, y1, z1, -lv.x, -lv.y, -lv.z, finlit, ld,
		 factor, FALSE,
		 1.0, 0.0, 0.0, 0.0,
		 0, &lred, &lgrn, &lblu);
	tmp_dred += cosln * lred;
	tmp_dgrn += cosln * lgrn;
	tmp_dblu += cosln * lblu;
      }

      if (hilit_sw) {
	double hx, hy, hz;

	hx = vx + lv.x;
	hy = vy + lv.y;
	hz = vz + lv.z;
	normalize(w, hx, hy, hz);
	coshn = -hx * nx - hy * ny - hz * nz;
	cosshn = sq(coshn);
	w = (cosshn - 1.0) / (cosshn * sq(attr->hilit_diff));
	if (w >= min_expseed) {
	  hl = exp(w);
	  coshv = hx * vx + hy * vy + hz * vz;
	  w = 2.0 * coshn * (cosln <= cosvn ? cosln : cosvn);
	  if (w < coshv) hl *= w / coshv;
	  if (hl > 0.0) {
	    tmp_sred += hl * lred;
	    tmp_sgrn += hl * lgrn;
	    tmp_sblu += hl * lblu;
	  }
	}
      }
    }

    ray_red += (cosvn * attr->lumi_red
		+ tmp_dred * diff_red
		+ tmp_sred * attr->hilit * attr->spec_red)
      * en_red;

    ray_grn += (cosvn * attr->lumi_grn
		+ tmp_dgrn * diff_grn
		+ tmp_sgrn * attr->hilit * attr->spec_grn)
      * en_grn;

    ray_blu += (cosvn * attr->lumi_blu
		+ tmp_dblu * diff_blu
		+ tmp_sblu * attr->hilit * attr->spec_blu)
      * en_blu;
  }

  if (++level <= trace_count) {
    double ks;
    double next_er, next_eg, next_eb;

    if (attr->tranflag) {   /* transparency object */
      double n2, nr;
      double px, py, pz;   /* refraction vector of V */
      double ar2, ag2, ab2;

      if (enter) {
	if (obj_complex_flag) calc_refract(object, x1f, y1f, z1f, &n2, &ar2, &ag2, &ab2);
	else {
	  n2 = attr->refract;
	  ar2 = attr->tran_red;
	  ag2 = attr->tran_grn;
	  ab2 = attr->tran_blu;
	}
      }
      else {
	CLUSTER *out_object;
	ANDFACT *out_factor;

	if (whereis(x1, y1, z1, vx, vy, vz, factor, &out_object, &out_factor)) {
	  if (obj_complex_flag)
	    calc_refract(out_object, x1f, y1f, z1f, &n2, &ar2, &ag2, &ab2);
	  else {
	    ATTR *at;

	    n2 = (at = out_factor->attribute)->refract;
	    ar2 = at->tran_red;
	    ag2 = at->tran_grn;
	    ab2 = at->tran_blu;
	  }
	}
	else {
	  n2 = 1.0;
	  ar2 = ag2 = ab2 = 0.0;
	}
      }
      nr = n2 / n1;

      if ((w = sq(nr) + sq(cosvn) - 1.0) >= 0.0) {
	double ncosvn, cosmnp, ncosmnp, ss, sp, kt;

	w = sqrt(w) - cosvn;
	px = (vx - nx * w) / nr;
	py = (vy - ny * w) / nr;
	pz = (vz - nz * w) / nr;

	ncosvn = nr * cosvn;
	ncosmnp = nr * (cosmnp = -nx * px - ny * py - nz * pz);
	ss = (cosvn - ncosmnp) / (cosvn + ncosmnp);
	sp = (ncosvn - cosmnp) / (ncosvn + cosmnp);
	kt = 1.0 - (ks = (sq(ss) + sq(sp)) * 0.5);

	if ((next_er = en_red * kt) >= min_energy ||
	    (next_eg = en_grn * kt) >= min_energy ||
	    (next_eb = en_blu * kt) >= min_energy)
	  trace(x1, y1, z1, px, py, pz, factor, enter, level, next_er, next_eg, next_eb, n2, ar2, ag2, ab2);
      }
      else ks = 1.0;
    }
    else ks = 1.0;

    next_er = en_red * ks * attr->spec_red;
    next_eg = en_grn * ks * attr->spec_grn;
    next_eb = en_blu * ks * attr->spec_blu;
    if (next_er >= min_energy || next_eg >= min_energy || next_eb >= min_energy) {
      w = 2.0 * cosvn;
      trace(x1, y1, z1, vx + nx * w, vy + ny * w, vz + nz * w, factor, enter, level, next_er, next_eg, next_eb, n1, ar1, ag1, ab1);
    }
  }

  return 255;
}

int sample(x, y)
double x, y;
{
  double w, vx, vy, vz;
  int alpha;

  vx = ca11 * x + ca21 * y + ca31z;
  vy = ca12 * x + ca22 * y + ca32z;
  vz = ca13 * x + ca23 * y + ca33z;
  normalize(w, vx, vy, vz);

  ray_red = ray_grn = ray_blu = 0.0;
  alpha = trace(camera_x, camera_y, camera_z, vx, vy, vz,
		NULL, FALSE,
		1,
		1.0, 1.0, 1.0,
		cam_n, cam_ar, cam_ag, cam_ab
		);
  sample_red = ray_red <= 0.0 ? 0 : (ray_red >= 1.0 ? 255 : ray_red * 256.0);
  sample_grn = ray_grn <= 0.0 ? 0 : (ray_grn >= 1.0 ? 255 : ray_grn * 256.0);
  sample_blu = ray_blu <= 0.0 ? 0 : (ray_blu >= 1.0 ? 255 : ray_blu * 256.0);

  return alpha;
}

/*
	compute the shade value for each pixel of a scan line
*/

void compute()
{
  double spp = sq((double)super_level);
  int px, py;
  double pitch_x, pitch_y;
  double x0, y0;


  pitch_x = pixel_size_x / super_level;
  pitch_y = pixel_size_y / super_level;

  y0 = screen_upper_edge - pitch_y / 2;

  for (py = 0; py < reso_v; py++) {
    char *p = scan_line_buffer;

    x0 = screen_left_edge + pitch_x / 2;

    for (px = 0; px < reso_h; px++) {
      int sx, sy;
      double x, y;
      double a, r, g, b;

      ray_red = ray_grn = ray_blu = 0.0;

      y = y0;

      for (sy = 0; sy < super_level; sy++) {

	x = x0;

	for (sx = 0; sx < super_level; sx++) {
	  double w, vx, vy, vz;
	  vector v;

	  vx = ca11 * x + ca21 * y + ca31z;
	  vy = ca12 * x + ca22 * y + ca32z;
	  vz = ca13 * x + ca23 * y + ca33z;
	  normalize(w, vx, vy, vz);
	  (void)trace(camera_x, camera_y, camera_z, vx, vy, vz,
		      NULL, FALSE,
		      1,
		      1.0, 1.0, 1.0,
		      cam_n, cam_ar, cam_ag, cam_ab
		      );
	  x += pitch_x;
	}
	y += pitch_y;
      }

      *p++ = 0;
      *p++ = (ray_red /= spp) <= 0.0 ? 0 : (ray_red >= 1.0 ? 255 : ray_red * 256.0);
      *p++ = (ray_grn /= spp) <= 0.0 ? 0 : (ray_grn >= 1.0 ? 255 : ray_grn * 256.0);
      *p++ = (ray_blu /= spp) <= 0.0 ? 0 : (ray_blu >= 1.0 ? 255 : ray_blu * 256.0);

      x0 += pixel_size_x;
    }
    y0 -= pixel_size_y;

    fwrite(scan_line_buffer, 4, reso_h, wfp);
  }
}

bool is_edge(a11, r11, g11, b11, a13, r13, g13, b13,
	     a31, r31, g31, b31, a33, r33, g33, b33)
int a11, a13, a31, a33,
    r11, g11, b11, r13,
    g13, b13, r31, g31,
    b31, r33, g33, b33;
{
  if (a11 && a13 && a31 && a33) {
    if ((r11 > r13 ? r11 - r13 : r13 - r11) > anti_thres ||
	(g11 > g13 ? g11 - g13 : g13 - g11) > anti_thres ||
	(b11 > b13 ? b11 - b13 : b13 - b11) > anti_thres)
      return TRUE;

    if ((r11 > r31 ? r11 - r31 : r31 - r11) > anti_thres ||
	(g11 > g31 ? g11 - g31 : g31 - g11) > anti_thres ||
	(b11 > b31 ? b11 - b31 : b31 - b11) > anti_thres)
      return TRUE;

    if ((r11 > r33 ? r11 - r33 : r33 - r11) > anti_thres ||
	(g11 > g33 ? g11 - g33 : g33 - g11) > anti_thres ||
	(b11 > b33 ? b11 - b33 : b33 - b11) > anti_thres)
      return TRUE;

    if ((r13 > r31 ? r13 - r31 : r31 - r13) > anti_thres ||
	(g13 > g31 ? g13 - g31 : g31 - g13) > anti_thres ||
	(b13 > b31 ? b13 - b31 : b31 - b13) > anti_thres)
      return TRUE;

    if ((r13 > r33 ? r13 - r33 : r33 - r13) > anti_thres ||
	(g13 > g33 ? g13 - g33 : g33 - g13) > anti_thres ||
	(b13 > b33 ? b13 - b33 : b33 - b13) > anti_thres)
      return TRUE;

    return
      ((r31 > r33 ? r31 - r33 : r33 - r31) > anti_thres ||
       (g31 > g33 ? g31 - g33 : g33 - g31) > anti_thres ||
       (b31 > b33 ? b31 - b33 : b33 - b31) > anti_thres);
  }
  else return (a11 || a13 || a31 || a33);
}

void anti_alias(x, y, depth, pitch_x, pitch_y,
		a11, r11, g11, b11, a13, r13, g13, b13,
		a31, r31, g31, b31, a33, r33, g33, b33,
		a, r, g, b)
double x, y;
int depth;
double pitch_x, pitch_y;
int a11, a13, a31, a33,
    r11, r13, r31, r33,
    g11, g13, g31, g33,
    b11, b13, b31, b33;
int *a, *r, *g, *b;
{
  if (depth < anti_depth &&
      is_edge(a11, r11, g11, b11, a13, r13, g13, b13,
	      a31, r31, g31, b31, a33, r33, g33, b33))
    {
      double pitch_xh, pitch_yh, x2, y2;

      int
	a12, a21, a22, a23, a32,
	r12, r21, r22, r23, r32,
	g12, g21, g22, g23, g32,
	b12, b21, b22, b23, b32;

      int
	pa1, pa2, pa3, pa4,
	pr1, pr2, pr3, pr4,
	pg1, pg2, pg3, pg4,
	pb1, pb2, pb3, pb4;

      x2 = x + (pitch_xh = pitch_x * 0.5);
      y2 = y + (pitch_yh = pitch_y * 0.5);
      a12 = sample(x2, y);
      r12 = sample_red;
      g12 = sample_grn;
      b12 = sample_blu;
      a21 = sample(x, y2);
      r21 = sample_red;
      g21 = sample_grn;
      b21 = sample_blu;
      a22 = sample(x2, y2);
      r22 = sample_red;
      g22 = sample_grn;
      b22 = sample_blu;
      a23 = sample(x + pitch_x, y2);
      r23 = sample_red;
      g23 = sample_grn;
      b23 = sample_blu;
      a32 = sample(x2, y + pitch_y);
      r32 = sample_red;
      g32 = sample_grn;
      b32 = sample_blu;
      anti_alias(x, y, ++depth, pitch_xh, pitch_yh,
		 a11, r11, g11, b11, a12, r12, g12, b12,
		 a21, r21, g21, b21, a22, r22, g22, b22,
		 &pa1, &pr1, &pg1, &pb1);
      anti_alias(x2, y, depth, pitch_xh, pitch_yh,
		 a12, r12, g12, b12, a13, r13, g13, b13,
		 a22, r22, g22, b22, a23, r23, g23, b23,
		 &pa2, &pr2, &pg2, &pb2);
      anti_alias(x, y2, depth, pitch_xh, pitch_yh,
		 a21, r21, g21, b21, a22, r22, g22, b22,
		 a31, r31, g31, b31, a32, r32, g32, b32,
		 &pa3, &pr3, &pg3, &pb3);
      anti_alias(x2, y2, depth, pitch_xh, pitch_yh,
		 a22, r22, g22, b22, a23, r23, g23, b23,
		 a32, r32, g32, b32, a33, r33, g33, b33,
		 &pa4, &pr4, &pg4, &pb4);
      *a = (pa1 + pa2 + pa3 + pa4) / 4;
      *r = (pr1 + pr2 + pr3 + pr4) / 4;
      *g = (pg1 + pg2 + pg3 + pg4) / 4;
      *b = (pb1 + pb2 + pb3 + pb4) / 4;
    }
  else
    {
      *a = a11 == 0 ? 0.0 : 0.25;
      *a += a13 == 0 ? 0.0 : 0.25;
      *a += a31 == 0 ? 0.0 : 0.25;
      *a += a33 == 0 ? 0.0 : 0.25;
      *r = (r11 + r13 + r31 + r33) / 4;
      *g = (g11 + g13 + g31 + g33) / 4;
      *b = (b11 + b13 + b31 + b33) / 4;
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
  static int mb[8] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};

  ptr += 7;

  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) putchar((*ptr & mb[j]) == 0 ? '0' : '1');
    --ptr;
  }

  putchar('\n');
}

#endif

void cam_refract()
{
  CLUSTER *out_object;
  ANDFACT *out_factor;

  if (cam_in_obj_flag) {
    double vx, vy, vz, w;

    w = sfabs(focal_dist);
    vx = ca31z / w;
    vy = ca32z / w;
    vz = ca33z / w;

    if (whereis(camera_x, camera_y, camera_z, vx, vy, vz, NULL, &out_object, &out_factor)) {
      if (obj_complex_flag)
	calc_refract(out_object, camera_x + vx * MARGIN,
		     camera_y + vy * MARGIN,
		     camera_z + vz * MARGIN,
		     &cam_n, &cam_ar, &cam_ag, &cam_ab);
      else {
	ATTR *at;

	cam_n = (at = out_factor->attribute)->refract;
	cam_ar = at->tran_red;
	cam_ag = at->tran_grn;
	cam_ab = at->tran_blu;
      }
    }
    else {
      cam_n = 1.0;
      cam_ar = cam_ag = cam_ab = 0.0;
    }
  }
  else {
    cam_n = 1.0;
    cam_ar = cam_ag = cam_ab = 0.0;
  }
}
