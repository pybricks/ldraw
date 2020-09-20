/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague and Don Heyse
 */

#include <math.h>
#include <float.h>

#include "lsynthcp.h"
#include "band.h"
#include "hose.h"

/*
 * 0 SYNTH BEGIN DEFINE BAND <fill> RUBBER_BAND "Descr" <scale> <thresh>
 * 1 <len>  a b c  d e f  g h i  j k l "name"
 * 1 <len>  a b c  d e f  g h i  j k l "name"
 * 0 SYNTH END
 */

int           n_band_types;
band_attrib_t   band_types[32];

#define N_BAND_TYPES n_band_types

/*
 * 0 SYNTH BEGIN DEFINE BAND CONSTRAINTS
 * 1 <dia>  a b c  d e f  g h i  j k l  "name"
 * 0 SYNTH END
 */

int    n_band_constraints;
part_t   band_constraints[64];

#define N_BAND_CONSTRAINTS n_band_constraints

// Make current band_type a global var so we can use it everywhere.
band_attrib_t *band_type = NULL;

/************************************************************************/
// Return 1 if the v2 bends left of v1, -1 if right, 0 if straight ahead.
int turn_vector(PRECISION v1[3], PRECISION v2[3])
{
  /* Pos for left bend, 0 = linear */
  PRECISION vec_product = (v1[0] * v2[1]) - (v1[1] * v2[0]);

  if (vec_product > 0.0) return(1);
  if (vec_product < 0.0) return(-1);
  return(0);
}

//**********************************************************************
void band_ini(void)
{
  int i;

  for (i = 0; i < N_BAND_TYPES; i++) {
    printf("%-20s = SYNTH BEGIN %s 16\n",band_types[i].type, band_types[i].type);
  }
}

void
list_band_types(void)
{
  int i;

  printf("\n\nBand type synthesizable parts\n");
  for (i = 0; i < N_BAND_TYPES; i++) {
    printf("  %-20s %s\n",band_types[i].type, band_types[i].descr);
  }
}

int
isbandtype(char *type)
{
  int i;

  for (i = 0; i < N_BAND_TYPES; i++) {
    if (strncasecmp(band_types[i].type,type,strlen(band_types[i].type)) == 0) {
      return 1;
    }
  }
  return 0;
}

int
isbandconstraint(char *type)
{
  int i;

  for (i = 0; i < N_BAND_CONSTRAINTS; i++) {
    if (strcasecmp(band_constraints[i].type,type) == 0) {
      return 1;
    }
  }
  return 0;
}

void
list_band_constraints(void)
{
  int i;

  printf("\n\nBand type synthesis constraints\n");
  for (i = 0; i < N_BAND_CONSTRAINTS; i++) {
    printf("    %11s\n",band_constraints[i].type);
  }
}

/*
 * Calculate the issues of crossers
 */

void
calc_crosses(
  LSL_band_constraint *k,
  LSL_band_constraint *m,
  int                 *layer,
  FILE *output)
{
  PRECISION xlk = k->end_line[0] - k->start_line[0];
  PRECISION ylk = k->end_line[1] - k->start_line[1];
  PRECISION xnm = m->end_line[0] - m->start_line[0];
  PRECISION ynm = m->end_line[1] - m->start_line[1];
  PRECISION xmk = m->start_line[0] - k->start_line[0];
  PRECISION ymk = m->start_line[1] - k->start_line[1];

  PRECISION det = xnm*ylk - ynm*xlk;
  if (fabs(det) < ACCY) {
    /* parallel lines */
  } else {
    PRECISION detinv = 1.0/det;
    PRECISION s = (xnm*ymk - ynm*xmk)*detinv;
    PRECISION t = (xlk*ymk - ylk*xmk)*detinv;
    if (s >= 0 && s <= 1.0 && t >= 0 && t <= 1.0) {
      PRECISION x = k->start_line[0] + xlk*s;
      PRECISION y = k->start_line[1] + ylk*s;
      if (k->n_crossings < 8) {
        k->crossings[k->n_crossings][0] = x;
        k->crossings[k->n_crossings][1] = y;
        k->crossings[k->n_crossings][2] = 0;
        k->n_crossings++;
      }
      if (m->n_crossings < 8) {
        m->crossings[m->n_crossings][0] = x;
        m->crossings[m->n_crossings][1] = y;
        m->crossings[m->n_crossings][2] = 0;
        m->n_crossings++;
      }
      if (k->layer == -1) {
        k->layer = *layer;
        *layer = *layer + 1;
      }
      if (m->layer == -1) {
        m->layer = *layer;
        *layer = *layer + 1;
      }
    }
  }
}

/*
 * Calculate the entry and exit angles of a band
 * around a constraint.
 */

void
calc_angles(
  band_attrib_t *type,
  LSL_band_constraint *k,
  FILE *output)
{
  PRECISION first_x, first_y, last_x, last_y;
  PRECISION dx,dy;
  PRECISION angle,ta;
  int i;
  float n;
  PRECISION pi = 2*atan2(1,0);

  if (k->cross || ! k->inside) {
    first_x = k->end_angle[0];
    first_y = k->end_angle[1];
    last_x  = k->start_angle[0];
    last_y  = k->start_angle[1];
  } else {
    first_x = k->start_angle[0];
    first_y = k->start_angle[1];
    last_x  = k->end_angle[0];
    last_y  = k->end_angle[1];
  }

  dx = (first_x - k->part.offset[0]);
  dy = (first_y - k->part.offset[1]);
  dx /= k->radius;
  dy /= k->radius;

  // Warning!!  acos() will give NAN if we give it badly normalized numbers.
  if (dx > 1.0) 
    dx = 1.0;
  if (dx < -1.0) 
    dx = -1.0;

  if (dy > 0) {
    angle = acos(dx);
  } else {
    angle = -acos(dx);
  }

#define REORIENT_TREAD 1
//#define SHOW_XY_PLANE_FOR_DEBUG 1
#define DEBUGGING_FIXED3_BANDS 1
// #define STRETCH_FIXED3 1

#define USE_TURN_ANGLE 1
#ifdef USE_TURN_ANGLE 
  // The problem here is this:
  // angle is the angle at which the tangent leaving this constraint starts at.
  // I need to get THE DIFFERENCE between that angle and the previous constraint
  // before I can calculate n_steps.
  // See get_turn_mat() fn in curve.t for calculating the turn angle

 {
  PRECISION r;
  PRECISION a[3];
  PRECISION b[3];

  extern PRECISION dotprod(PRECISION a[3], PRECISION b[3]);
  
  if (k->cross || ! k->inside) {
    vectorcp(a, k->end_angle);
    vectorcp(b, k->start_angle);
  } else {
    vectorcp(b, k->end_angle);
    vectorcp(a, k->start_angle);
  }

  vectorsub(a, k->part.offset);
  vectorsub(b, k->part.offset);
  normalize(a);
  normalize(b);

    //Dot product gives turn angle.  a.b=|a||b|cos(theta)
    //We normalized so |a|=|b|=1, which means theta = acos(a.b).
    r = dotprod(b, a);
    // Warning!!  acos() will give NAN if we give it badly normalized numbers.
    if (r > 1.0) 
      r = 1.0;
    if (r < -1.0) 
      r = -1.0;

    ta = r;
    r = acos(r);

    // Check crossprod(b, a) to see if r is CW or CCW.
    i = turn_vector(b, a); 

    // Handle round off error near 0 degree turn
    if ((ta + ACCY) > 1.0) // dotprod() is really 1 so call it no turn.
      i = 0;

    if (i > 0) 
      r = 2*pi - r;    

    else if (i == 0) // Handle 0 or 180 by the sign of dotprod()
    {
      if (ta < 0)
	r = pi;
      else if (k->layer == -2) // The ONLY constraint.
	r = 2*pi; // Cover the whole thing.
      else
	r = 0; // Skip on past this constraint.
    }

    if (k->layer == -2)
      k->layer == -1;

#ifdef DEBUGGING_FIXED3_BANDS
  if (k->cross || ! k->inside) 
    printf("OUT(%.2fx, %.2fy, %dr)  A = %.2f from (%.2f, %.2f) r = %.2f (%dT%.2f)\n", 
	   k->part.offset[0], k->part.offset[1], (int)k->radius,
	   angle * 180 / pi, dx, dy, r * 180 / pi, i, ta);
  else
    printf("IN (%.2fx, %.2fy, %dr)  A = %.2f from (%.2f, %.2f) r = %.2f (%dT%.2f)\n", 
	   k->part.offset[0], k->part.offset[1], (int)k->radius, 
	   angle * 180 / pi, dx, dy, r * 180 / pi, i, ta);
#endif

  ta = r;
 }

#endif

  k->s_angle = angle;

  // NOTE: Start converting FIXED3 to FIXED(N) with N segments.
  // Probably should use something more like a convex hull for FIXED(N)
  // Gotta review all this STRETCH_FIXED3 stuff.  What was I thinking???
  //
  // I should measure the length of the convex hull and stretch it to 
  // fit N segments if needed.  (Maybe also compress it if a bit too big?)
  //
  // Remember FIXED3 was created to indicate we need 3 part types to describe
  // a (probably rubbery) band (straights, arcs, and transitions).  
  // That still applies, but now we also want to know how many segments,
  // so extend it to FIXED(N) like the FIXED(N) hoses.
  // After all, when you think about it, you only need the special blending
  // parts when the material is flexible but segmented, which almost implies  
  // rubber treads.  However I suppose we could retain the special FIXED3
  // function just in case they make something in multiple sizes.
  // We don't actually lose much even if we store FIXED3 as type->fill == 3
  // because a 3 sided rubber tread is pretty dull stuff.
  // But storing it that way does simplify a bunch of code.
  // 
  // Now, do we need full convex hull code or can we get by with some 
  // turn_vector() tests.  We can borrow that fn from stub.c in ldglite.
  // Are we always going CCW around the constraints?

  //printf("Band FillType = %d \n",type->fill);
  // If (type->fill > FIXED) it contains the length (in segments) of the band.

  if ((type->fill == FIXED3) || (type->fill > FIXED)) {
    PRECISION circ;
    if (angle < 0) {
      angle = - angle;
    }
#ifdef USE_TURN_ANGLE 
    circ = ta*k->radius;
#else
    circ = angle*2*k->radius;
#endif
#ifdef STRETCH_FIXED3
    // Do not round up the number of arc steps.  Stretch the tangent lines instead.
    k->n_steps = circ*type->scale;
#else
    k->n_steps = circ*type->scale+0.5;
#endif
  } else if (type->fill == FIXED) {

#ifdef USE_TURN_ANGLE 
    PRECISION circ;
    circ = ta*k->radius;
    k->n_steps = circ*type->scale + 0.5;
    k->n_steps++; // Not really steps, but segment endpoints?  So add 1 more point.
    printf("nsteps = %d = (%.2f / %.2f\n", k->n_steps, circ, 1.0 / type->scale);
#else
    n = type->scale * 2 * pi * k->radius + 0.5;

    // circumference
    for (i = 0; i < n; i++) {
      PRECISION f;
      f = i;
      f /= n;
      ta = angle - 2*pi*(1-f);
      dx = k->radius*cos(ta) + k->part.offset[0] - last_x;
      dy = k->radius*sin(ta) + k->part.offset[1] - last_y;
      if (sqrt(dx*dx+dy*dy) < type->thresh) {
        break;
      }
    }
    k->n_steps = i+1;
#endif
  } else { // (type->fill == STRETCH)

    n =  2 * pi * k->radius/band_res + 0.5;

    // circumference
    for (i = 0; i < n; i++) {
      PRECISION f;
      f = i;
      f /= n;
      ta = angle - 2*pi*(1-f);
      dx = k->radius*cos(ta) + k->part.offset[0] - last_x;
      dy = k->radius*sin(ta) + k->part.offset[1] - last_y;
      if (sqrt(dx*dx+dy*dy) < type->thresh) {
        break;
      }
    }
    k->n_steps = i+1;

  }
}

/*
 * figure out the intersections of a line and a circle.
 *   - line is in x = xo + f*t, y = yo + g*t (normalized parametric) form 
 *   - circle is (xj, yj, rj)
 */

int intersect_line_circle_2D(
  PRECISION xo, 
  PRECISION yo,
  PRECISION f,
  PRECISION g,
  PRECISION xj,
  PRECISION yj,
  PRECISION rj,
  PRECISION *x,
  PRECISION *y)
{
  PRECISION fsq, gsq, fgsq;
  PRECISION xjo, yjo;
  PRECISION fygx;
  PRECISION fxgy;
  PRECISION root;
  PRECISION t;

  fsq = f * f;
  gsq = g * g;
  fgsq = fsq + gsq; // dx2 + dy2 (should be 1 if normalized).

  if (fgsq < ACCY) {
    printf("line coefficients are corrupt\n"); // Because dx = dy = 0.
  }

  xjo = xj - xo;
  yjo = yj - yo;
  fygx = f*yjo - g*xjo;
  root = rj*rj*fgsq - fygx*fygx;

  // What's the story with this message?  I get it sometimes but this still works.
  // We eliminate the two degenerate cases (circles overlap) in calc_tangent_line() 
  // So this fn should ALWAYS work when it's called.  I suspect this is bogus...
  if (root < -ACCY) {
    printf("line does not intersect with circle\n");
  }

  fxgy = f*xjo + g*yjo;

  t = fxgy/fgsq;

  *x = xo + f*t;
  *y = yo + g*t;
  return 0;
}

/*
 * determine the tangent we want.
 */

int calc_tangent_line(
  LSL_band_constraint *k,
  LSL_band_constraint *l,
  FILE                *output)
{
  int inside1, inside2;
  PRECISION rl,rk,rlk;
  PRECISION xlk,xlksq,ylk,ylksq;
  PRECISION denom;
  PRECISION radius;
  PRECISION angle;
  PRECISION rx,ry;

  inside1 = k->inside;
  inside2 = l->inside;

  if (l->was_cross) {
    if (l->cross) {
      inside2 ^= 1;
    } else {
      inside1 ^= 1;
    }
  }

  rl = l->radius;
  rk = k->radius;

  switch ((inside1 << 1) | inside2) {
    case 3: /* inside to inside */
      /* end angle for i , start for j */
      /* start point line of ij line stored in j */
      /* endpoint of line of ij line stored in j */
    break;
    case 2: /* inside to outside */
      rl = - rl;
    break;
    case 1: /* outside to inside */
      rk = -rk;
    break;
    case 0: /* outside to outside */
      rk = -rk;
      rl = -rl;
    break;
  }

  rlk = rl - rk;

  xlk = l->part.offset[0] - k->part.offset[0];
  ylk = l->part.offset[1] - k->part.offset[1];

  xlksq = xlk*xlk;
  ylksq = ylk*ylk;

  denom = xlksq + ylksq;

  if (denom < ACCY) {
    /* circles are coincident - badness */
  } else {
    PRECISION root;

    root = denom - rlk*rlk;
    if (root < -ACCY) {
      /* tangent doesn exist */
      // ie. sqr(dist) < sqr(radius)
      // NOTE:  actually we should check for (root < +ACCY) because:
      // For the outer tangent case this means one circle is completely
      // inscribed within the other.
      // For the inner tangent case this means the circles touch or
      // overlap.
    } else {

      PRECISION a,b,c;
      PRECISION deninv,factor;
      PRECISION xo,yo;
      PRECISION f,g;

      if (root < 0) {
        root = 0;
      }
      root = sqrt(root);
      deninv = 1.0/denom;
      a = (-rlk*xlk - ylk*root)*deninv;
      b = (-rlk*ylk + xlk*root)*deninv;
      c = -(rk + a*k->part.offset[0] + b*k->part.offset[1]);

      /* we have the normalized form of the tangent line */

      /* now we map the line to parametric form */
      root = 1.0/(a * a + b * b);
      factor = -c*root;
      xo = a*factor;
      yo = b*factor;

      root = sqrt(root);

      f =  b*root;
      g = -a*root;

      /* now line is in x = xo + f*t, y = yo + g*t form */
      /* calculate endpoints of each line */

      intersect_line_circle_2D(
        xo,
        yo,
        f,
        g,
        k->part.offset[0],
        k->part.offset[1],
        k->radius,
       &k->start_line[0],
       &k->start_line[1]);
        k->start_line[2] = 0;
      vectorcp(k->start_angle,k->start_line);

      intersect_line_circle_2D(
        xo,
        yo,
        f,
        g,
        l->part.offset[0],
        l->part.offset[1],
        l->radius,
       &k->end_line[0],
       &k->end_line[1]);
        k->end_line[2] = 0;
      vectorcp(l->end_angle,k->end_line);

      // this means we need our previous neighbor's end line and our
      // start line to know the arc
    }
  }
  return 0;
}

/*
 * Render the arc around the constraint and the line to the next
 * constraint
 */

int draw_arc_line(
  band_attrib_t       *type,
  LSL_band_constraint *constraint,
  int                  color,
  int                  draw_line,
  FILE                *output,
  int                  ghost,
  char                *group,
  part_t      *absolute,
  LSL_band_constraint *f_constraint)
{
  int       i,j,n;
  PRECISION dx,dy,dz;
  PRECISION L1,L2;
  PRECISION inv[3][3];
  int steps;

  if (draw_line) {

    for (j = 0; j < constraint->n_crossings - 1; j++) {
      PRECISION orient[3][3];

      // determine the orientation of the part in the XY plane

      dx = constraint->crossings[j+1][0] - constraint->crossings[j][0];
      dy = constraint->crossings[j+1][1] - constraint->crossings[j][1];
      dz = constraint->crossings[j+1][2] - constraint->crossings[j][2];

      L1 = sqrt(dx*dx + dy*dy);
      L2 = sqrt(dx*dx + dy*dy + dz*dz);
#ifdef DEBUGGING_FIXED3_BANDS
      printf("Direction = (%.3f, %.3f, %.3f) => %.3f_L1, %.3f_L2)\n", dx, dy, dz, L1, L2);
      // ******************************
      // Based on this it looks like L2 is already in whole number units of part len
      // (where part len is 1/type->scale).
      // In order to be able to stretch things a bit to fit, I need a raw length
      // in L1 or L2.  (L2 should always = L1 since dz=0 when we move to XY plane.)
      //
      // Maybe not.  It just works out that way here in my example.
      // Really what I need to do is delay the line drawing until after the arcs.
      // Calculate the arcs (making sure to always round n arc segs down)
      // Then move the j+1 crossings to meet the ends of the arcs, 
      // stretching L1 a bit in the process.
#endif
      if (L1 == 0) {

        orient[0][0] = 1;
        orient[1][0] = 0;
        orient[2][0] = 0;
        orient[0][1] = 0;
        orient[1][1] = 0;
        orient[2][1] =-1;
        orient[0][2] = 0;
        orient[1][2] = 1;
        orient[2][2] = 0;
      } else {
        orient[0][0] =  dy/L1;  //  cos
        orient[1][0] = -dx/L1;  //  sin
        orient[2][0] =  0;
        orient[0][1] =  dx/L2;  // -sin
        orient[1][1] =  dy/L2;  //  cos
        orient[2][1] =  dz/L2;
        orient[0][2] = -dx*dz/(L1*L2);
        orient[1][2] = -dy*dz/(L1*L2);
        orient[2][2] =  L1/L2;
      }

      if (type->fill == STRETCH) {
        n = 1;
        steps = 0;
      } else if (type->fill == FIXED3) {
        n = L2*type->scale+0.5;
        steps = 1;
      } else if (type->fill > FIXED) {
        n = L2*type->scale+0.5;
        steps = 1;
      } else { // FIXED
        n = L2*type->scale+0.5;
        steps = 0;
      }
      
#ifdef STRETCH_FIXED3
      L1 = 1;
#endif
      for (i = steps; i < n; i++) {
        part_t part;
        PRECISION tm[3][3];
        PRECISION foffset[3];

        part = type->tangent;

        if (type->fill == STRETCH) {
          PRECISION scale[3][3];
          int i,j;

          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              scale[i][j] = 0;
            }
          }
          scale[0][0] = 1;
          scale[1][1] = L2;
          scale[2][2] = 1;
          matrixmult(part.orient,scale);
        }
#ifdef STRETCH_FIXED3
	else if ((type->fill == FIXED3) || (type->fill > FIXED)) {
          PRECISION scale[3][3];
          int i,j;

          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              scale[i][j] = 0;
            }
          }
          scale[0][0] = 1;
	  scale[1][1] = 1;
          scale[2][2] = 1;
	  if (n > 1)
	    scale[0][0] = L1 = (L2*type->scale -1)/(n-1);
          matrixmult(part.orient,scale);
        }
	
#endif
#ifdef DEBUGGING_FIXED3_BANDS
	if (i == steps)
	  printf("Scale(%.2f, %.3f, %d) = %.2f)\n", L2, type->scale, n, L1);
#endif

        // We performed this to start:
        //   1.  move the assembly so the first constraint is at the origin
        //   2.  rotate offset ations about the inverse of the first constraint's
        //       orientation bringing everything into the X/Y plane

        // Now we're putting together the segments and putting them back into place
        //   1.  multiply tangent part orientation times the orient
        //       so we know how to orient the tangent part in the X/Y plane.
        //   2.  rotate the tangent offsets into the X/Y plane.

        matrixcp(tm,part.orient);
        matrixmult3(part.orient,orient,tm);

        vectorrot3(part.offset,type->tangent.offset,part.orient);

        //   3.  Calculate the tangent part's orientation in the X/Y plane

        foffset[0] = constraint->crossings[j][0] + dx * i / n - part.offset[0];
        foffset[1] = constraint->crossings[j][1] + dy * i / n - part.offset[1];
        foffset[2] = constraint->crossings[j][2] + dz * i / n - part.offset[2];

#ifdef STRETCH_FIXED3
        foffset[0] += (L1-1) * dx * (i-1) / ((PRECISION)n -0.5);
        foffset[1] += (L1-1) * dy * (i-1) / ((PRECISION)n -0.5);
        foffset[2] += (L1-1) * dz * (i-1) / ((PRECISION)n -0.5);
#endif

        vectorcp(part.offset,foffset);

#ifdef REORIENT_TREAD
        //   4.  Reorient the tangent part in the X/Y plane based on the first
        //       constraint's offset (for things like technic turntable top,
        //       where the gear plane does not go through the origin.

        vectorsub( part.offset,band_constraints[f_constraint->band_constraint_n].offset);

        //   5.  Orient tangent part based on the orientation of the first
        //       constraint's orientation (for things like technic turntable
        //       who's gear plane is perpendicular to the plane of say the 24T
        //       gears.

	matrixinv(inv,band_constraints[f_constraint->band_constraint_n].orient);
        vectorrot( part.offset,inv);

        //   6.  Orient the tangent part offsetation back to the absolute 3D
        //       offsetation.

        vectorrot( part.offset,absolute->orient);

        //   7.  Now add the absolute orientation of the first constraint

        vectoradd( part.offset,absolute->offset);

        //   8.  Change the tangent part orientation based on the first
        //       constraint's orientation (again for things like the
        //       technic turntable top, where the gear plane is perpendicular
        //       to the standard 24T gear's gear plane

        matrixcp(tm,part.orient);
        matrixmult3(part.orient, inv, tm);

        //   9.  Now rotate the part back into its correct orientation in 3D
        //       space.

        matrixcp(tm,part.orient);
        matrixmult3(part.orient,absolute->orient,tm);
#endif // REORIENT_TREAD

        output_line(
          output,
          ghost,
          group,
          color,
          part.offset[0],part.offset[1],part.offset[2],
          part.orient[0][0],part.orient[0][1],part.orient[0][2],
          part.orient[1][0],part.orient[1][1],part.orient[1][2],
          part.orient[2][0],part.orient[2][1],part.orient[2][2],
          type->tangent.type);
      }
    }
  }

  // Create the arc

  {
    PRECISION pi = 2*atan2(1,0);
    PRECISION f[3];

    /* now for the arc */

    if (type->fill == STRETCH) {
      n = 2*pi*constraint->radius/band_res;
    } else {
      n = 2*pi*constraint->radius*type->scale;
    }

    // vector for the first arc part

    f[0] = constraint->radius * cos(constraint->s_angle);
    f[1] = constraint->radius * sin(constraint->s_angle);
    f[2] = 0;

    if ((type->fill == FIXED3) || (type->fill > FIXED)) {
      steps = constraint->n_steps + 2;
    } else {
      steps = constraint->n_steps;
    }

    for (i = 1; i < steps; i++) {
      PRECISION orient[3][3];
      PRECISION foffset[3];
      PRECISION tm[3][3];
      PRECISION s[3];
      part_t    part;

      // for FIXED3 (e.g. rubber tread), the first and last segments are
      // treated as transition pieces, otherwise just use arc parts

#ifdef NEW_BAND_CODE_BUT_STILL_NEEDS_WORK
      if (((type->fill == FIXED3) || (type->fill > FIXED))
	  && i == 1 && i+1 == steps) {
	part = type->tangent; // Only one part.  It probably should be straight.
	// Of course the right way would be to check the turn angle here,
	// but we don't save it (r) from calc_angles().
      } else
#endif
      if (((type->fill == FIXED3) || (type->fill > FIXED)) && i == 1) {
        part = type->start_trans;
      } if (((type->fill == FIXED3) || (type->fill > FIXED)) && i+1 == steps) {
        part = type->end_trans;
      } else {
        part = type->arc;
      }

      // rotate the arc part so it hugs the current constraint

      s[0] = constraint->radius * cos(constraint->s_angle + 2*pi*i/n);
      s[1] = constraint->radius * sin(constraint->s_angle + 2*pi*i/n);
      s[2] = 0;

      dx = s[0] - f[0];
      dy = s[1] - f[1];
      dz = s[2] - f[2];

      L1 = sqrt(dx*dx + dy*dy);
      L2 = sqrt(dx*dx + dy*dy + dz*dz);

      if (L1 == 0) {
        orient[0][0] =  1;
        orient[1][0] =  0;
        orient[2][0] =  0;
        orient[0][1] =  0;
        orient[1][1] =  0;
        orient[2][1] = -1;
        orient[0][2] =  0;
        orient[1][2] =  1;
        orient[2][2] =  0;
      } else {
        orient[0][0] =  dy/L1;
        orient[1][0] = -dx/L1;
        orient[2][0] =  0;
        orient[0][1] =  dx/L2;
        orient[1][1] =  dy/L2;
        orient[2][1] =  dz/L2;
        orient[0][2] = -dx*dz/(L1*L2);
        orient[1][2] = -dy*dz/(L1*L2);
        orient[2][2] =  L1/L2;
      }

      if (type->fill == STRETCH) {
        PRECISION scale[3][3];
        PRECISION angle = 2 * pi / n;
        PRECISION l = type->scale*sin(angle);
        int k,j;

        if (i + 1 == steps) {
          l *= 5;
        }

        for (k = 0; k < 3; k++) {
          for (j = 0; j < 3; j++) {
            scale[k][j] = 0;
          }
        }
        scale[0][0] = 1;
        scale[1][1] = L2+l;
        scale[2][2] = 1;
        matrixmult(part.orient,scale);
      }

      // Now we're putting together the segments and putting them back into place
      //   1.  multiply arc part orientation times the orient
      //       so we know how to orient the arc part in the X/Y plane.
      //   2.  rotate the arc offsets into the X/Y plane.

      matrixcp(tm,part.orient);
      matrixmult3(part.orient,orient,tm);

#ifdef NEW_BAND_CODE_BUT_STILL_NEEDS_WORK
      if (((type->fill == FIXED3) || (type->fill > FIXED))
	  && i == 1 && i+1 == steps) // if (part == type->tangent)
        vectorrot3(part.offset,type->tangent.offset,part.orient);
      else
#endif
      vectorrot3(part.offset,type->arc.offset,part.orient);

      //   3.  Calculate the arc part's offsetation in the X/Y plane

      foffset[0] = f[0] + constraint->part.offset[0] - part.offset[0];
      foffset[1] = f[1] + constraint->part.offset[1] - part.offset[1];
#ifdef WTF_IS_THIS
      foffset[2] = f[2] + constraint->part.offset[2] - part.offset[2];
#else
      // I thought the constraint Z coordinate was always zero in the X/Y plane.
      foffset[2] = f[2] + 0 - part.offset[2];
#endif

      vectorcp(part.offset,foffset);

#ifdef REORIENT_TREAD
      //   4.  Rotate the arc part in the X/Y plane based on the first
      //       constraint's offset (for things like technic turntable top,
      //       where the gear plane does not go through the origin.)

      vectorsub(part.offset,band_constraints[f_constraint->band_constraint_n].offset);

      //   5.  Orient arc part based on the orientation of the first
      //       constraint's orientation (for things like technic turntable
      //       who's gear plane is perpendicular to the plane of say the 24T
      //       gears.)

      matrixinv(inv,band_constraints[f_constraint->band_constraint_n].orient);
      vectorrot(part.offset,inv);

      //   6.  Orient the arc part offsetation back to the absolute 3D
      //       offsetation.

      vectorrot(part.offset,absolute->orient);

      //   7.  Now add the absolute offsetation of the first constraint

      vectoradd( part.offset,absolute->offset);

      //   8.  Change the arc part orientation based on the first
      //       constraint's orientation (again for things like the
      //       technic turntable top, where the gear plane is perpendicular
      //       to the standard 24T gear's gear plane)

      matrixcp(tm,part.orient);
      matrixmult3(part.orient, inv, tm);

      //   9.  Now rotate the part back into its correct orientation in 3D
      //       space.

      matrixcp(tm,part.orient);
      matrixmult3(part.orient,absolute->orient,tm);
#endif // REORIENT_TREAD

      output_line(
        output,
        ghost,
        group,
        color,
        part.offset[0],part.offset[1],part.offset[2],
        part.orient[0][0],part.orient[0][1],part.orient[0][2],
        part.orient[1][0],part.orient[1][1],part.orient[1][2],
        part.orient[2][0],part.orient[2][1],part.orient[2][2],
        part.type);

      vectorcp(f,s);
    }
  }

  return 0;
}

void
showconstraints(
  FILE                *output,
  LSL_band_constraint *constraints,
  int                  n_constraints,
  int                  color)
{
#if 1
  int i;

  for (i = 0; i < n_constraints; i++) {
    part_t *cp = &constraints[i].part;
    output_line(
      output,
      0,
      NULL,
      color,
      cp->offset[0],   cp->offset[1],   cp->offset[2],
      cp->orient[0][0],cp->orient[0][1],cp->orient[0][2],
      cp->orient[1][0],cp->orient[1][1],cp->orient[1][2],
      cp->orient[2][0],cp->orient[2][1],cp->orient[2][2],
      cp->type);
  }

  fflush(output);
#endif
}

static void
rotate_constraints(
  LSL_band_constraint *constraints,
  int                  n_constraints,
  PRECISION            m[3][3])
{
  int i;
  PRECISION t[3][3];

  for (i = 0; i < n_constraints-1; i++) {
    vectorrot(constraints[i].part.offset,m);
    matrixcp(t,constraints[i].part.orient);
    matrixmult3(constraints[i].part.orient,m,t);
  }
}

/*
 * This subroutine synthesizes planar rubber bands, chain, and treads.
 *
 * We do all the arc and tangent analysis in the X/Y plane.
 *
 * The synthesis plane is defined by the first constraint.  First we move the
 * first constaint to the origin. Then we calculate the inverse of the first
 * constraints orientation, and multiply all the constraints' offsetations by
 * the inverse of the first constraint's orientation.
 *
 * Some of the gears' are oriented in the X/Y plane, while other gears are
 * oriented in the X/Z plane.  The band_constraints array above, describes each of the
 * supported constraint types, and their orientation.  We multiply all the
 * constraints by the first constraint's band_constraints orientation.
 *
 * Also, in some cases, some of the constraint types described in band_constraints
 * need to be offset to get the place where the band should hit onto the
 * X/Y plane.
 */

//*****************************************************************
// NOTES:  
// 
// What's the meaning of the tangent directives? (INSIDE, OUTSIDE, CROSS)
// INSIDE, OUTSIDE refer to where we want the next wheel to be.
// It's either gonna be placed inside the band, or outside the band.
// If I move the constraints in the XY plane and draw the band CCW
// then the next wheel will be left of the band for INSIDE.
// The constraint will be to the right of the band for OUTSIDE.
// The default is INSIDE.  
// CROSS means OUTSIDE, then INSIDE.  It can be used to make an X.
// 
// CROSS really only should be allowed for string and rubber bands.
// It doesn't make much sense for rubber treads, although I suspect
// it could be used to put one wheel of three on the outside.
//
//*****************************************************************

int
synth_band(
  char *type,
  int n_constraints,
  LSL_band_constraint *constraints,
  int color,
  FILE *output,
  int ghost,
  char *group)
{
  int i;
  int cross = 0;
  int was_cross = 0;
  int inside = 1;
  int first,last;
  part_t absolute;
  PRECISION inv[3][3],trot[3][3];
  int layer = 0;

  // Make band_type a global var so we can use it everywhere.
  // band_attrib_t *band_type = NULL;

  /* Search for band type */
  band_type = NULL;
  for (i = 0; i < N_BAND_TYPES; i++) {
    if (strcasecmp(type,band_types[i].type) == 0) {
      band_type = &band_types[i];
      break;
    }
  }
  if (band_type == NULL) {
    return 0;
  }

  // Hmmm, maybe I should redirect all info messages to stderr
  // (so lsynth can be used as a filter).
  // Or maybe redirect to stderr ONLY if its a filter (if (output == stdout))

  if (n_constraints < 1) {
    fprintf(stderr, "No BAND constraints found.\n");
    return 0;
  }

  first = -1;

  for (i = 0; i < n_constraints; i++) {
    constraints[i].radius       = 0;
    constraints[i].inside       = inside;
    constraints[i].cross        = cross;
    constraints[i].was_cross    = was_cross;
    constraints[i].n_crossings  = 0;
    constraints[i].layer        = -1;
    was_cross = 0;
    if (strcasecmp(constraints[i].part.type,"INSIDE") == 0) {
      inside = 1;
    } else if (strcasecmp(constraints[i].part.type,"OUTSIDE") == 0) {
      inside = 0;
    } else if (strcasecmp(constraints[i].part.type,"CROSS") == 0) {
      inside ^= 1;

    } else {
      int k;

      // search the constraints table

      for (k = 0; k < N_BAND_CONSTRAINTS; k++) {
        if (strcasecmp(constraints[i].part.type,band_constraints[k].type) == 0) {
          constraints[i].band_constraint_n = k;

          constraints[i].radius = band_constraints[k].attrib;
          break;
        }
      }
    }
    if (first == -1 && constraints[i].radius) {
      first = i;
    }
  }

#ifdef USE_TURN_ANGLE 
  if (n_constraints == 1) {
    constraints[first].layer        = -2; // Mark this as the ONLY constraint.
  }
#endif

  /* create an extra constraint that represents the final state of the
   * first pulley */
  if (band_type->pulley == 0) {
    memcpy(&constraints[i],&constraints[first],sizeof(constraints[i]));
    n_constraints++;
  }
  constraints[i].inside    = inside;
  constraints[i].cross     = cross;
  constraints[i].was_cross = was_cross;

  /* record the first constraint in its original form */

  absolute = constraints[first].part;

  /* 1. move the first constraint to the origin */

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      vectorsub(constraints[i].part.offset,absolute.offset);
    }
  }

  // showconstraints(output,constraints,n_constraints,14);

  /* 2. bring the entire assembly into the part's natural orientation */

  matrixinv(inv,absolute.orient);

  rotate_constraints(constraints,n_constraints,inv);

  // showconstraints(output,constraints,n_constraints,4);

  /* 3. bring the assembly into the X/Y plane (necessary for first constraints
   *    who's gear plane is different that the default gear plane used by
   *    simple gears, like technic turntable).
   */

  rotate_constraints(constraints,n_constraints,
    band_constraints[constraints[first].band_constraint_n].orient);

  //showconstraints(output,constraints,n_constraints,15);

  /* 4. Now that the whole assembly is in AN the X/Y plane, move everything
   *    so the center of the axle of the first constraint is at the origin.
   *    With the current list of band constraints that means maybe adjusting 
   *    the Z Coordinates.
   */

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      vectoradd(constraints[i].part.offset, //was vectorsub()
        band_constraints[constraints[first].band_constraint_n].offset);
    }
  }

#if 0
  /* 5. Offset constraints that were not modeled at the orgin (or on an axis).
   *    
   *    At this point we can assume all constraint centers are in THE X/Y plane
   *    at Z=0.  However the part origins may need to be moved in X and/or Y 
   *    if for some reason the part was not modeled symmetrically about either 
   *    the X, Y, or Z axis.   Do not change the Z at this point.  
   *    And skip the first constraint because we already did it above (step 4).
   *    
   *    Do not undo this step in the final reorient code.  It's one way.
   */

  for (i = 1; i < n_constraints; i++) { // Start at constraint 1, not 0.
    if (constraints[i].radius) {
      vectoradd(constraints[i].part.offset,
        band_constraints[constraints[i].band_constraint_n].offset);
      constraints[i].part.offset[2] = 0; // Set Z to zero, just in case.
    }
  }
#endif

  //***************************************************************************
  //NOTE: This dies if I use only constraints  NOT listed as band constraints.
  //      By the time I get here constraints[0].part->type is a mangled string.
  //      Perhaps the whole constraint part is crap, or perhaps just the name.
  //***************************************************************************

  fflush(output);
#ifdef DEBUGGING_FIXED3_BANDS
  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      printf("constraint[%d] = (%d, %d, %d)\n", i, 
	     (int)constraints[i].part.offset[0], 
	     (int)constraints[i].part.offset[1], 
	     (int)constraints[i].part.offset[2] );
    }
  }
#endif
#ifdef SHOW_XY_PLANE_FOR_DEBUG
  showconstraints(output,constraints,n_constraints,3);
#endif

  /* figure out the tangents' intersections with circles */

  first = -1;
  last  = -1;

  for (i = 0; i < n_constraints - 1; ) {
    if (constraints[i].radius) {
      int j;
      if (first == -1) {
        first = i;
      }
      for (j = i+1; j < n_constraints; j++) {
        if (constraints[j].radius) {
#ifdef DEBUGGING_FIXED3_BANDS
	  printf("calc_tan(%d->%d)\n", i, j);
#endif
          calc_tangent_line(&constraints[i],&constraints[j],output);
          i = j;
          last = j;
          break;
        }
      }
      if (j == n_constraints+1) {
        i++;
      }
    } else {
      i++;
    }
  }
  vectorcp(constraints[first].end_angle,constraints[last].end_angle);

  /* calculate intersections between band straight line segments, so
   * we can make the line segments go around each other.
   */

  for (i = 0; i < n_constraints; i++) {
    if (constraints[i].radius) {
      vectorcp(constraints[i].crossings[0],constraints[i].start_line);
      constraints[i].n_crossings = 1;
    }
  }
  for (i = 0; i < n_constraints- 1; i++) {
    if (constraints[i].radius) {
      int j;
      for (j = i+1; j < n_constraints; j++) {
        if (constraints[j].radius) {
          calc_crosses(&constraints[i],&constraints[j],&layer,output);
        }
      }
    }
  }


#ifdef SHOW_XY_PLANE_FOR_DEBUG
  for (i = 0; i < n_constraints-1; i++) {
    part_t *cp = &constraints[i].part;
    LSL_band_constraint *k = &constraints[i];
    int color = 4;
    output_line(
      output,
      0,
      NULL,
      color,
      //cp->offset[0]+k->start_line[0],cp-> offset[1]+k->start_line[0], cp->offset[2]+k->start_line[2],
      k->start_line[0], k->start_line[1], k->start_line[2],
      cp->orient[0][0],cp->orient[0][1],cp->orient[0][2],
      cp->orient[1][0],cp->orient[1][1],cp->orient[1][2],
      cp->orient[2][0],cp->orient[2][1],cp->orient[2][2],
      "LS02.dat");
    color = 2;
    output_line(
      output,
      0,
      NULL,
      color,
      //cp->offset[0]+k->start_line[0],cp-> offset[1]+k->start_line[0], cp->offset[2]+k->start_line[2],
      k->end_line[0], k->end_line[1], k->end_line[2],
      cp->orient[0][0],cp->orient[0][1],cp->orient[0][2],
      cp->orient[1][0],cp->orient[1][1],cp->orient[1][2],
      cp->orient[2][0],cp->orient[2][1],cp->orient[2][2],
      "LS02.dat");
  }
#endif


#define BAND_DIAM 4

  /* calculate the depth for each band at crossing */

  if (layer > 0) {
    {
      PRECISION layer_offset;
      layer_offset = (((layer) / 2) % 2) * BAND_DIAM/2;
      //layer_offset = 0;
      for (i = 0; i < n_constraints-1; i++) {
        if (constraints[i].radius) {
          int j;
          for (j = 1; j < constraints[i].n_crossings; j++) {
            constraints[i].crossings[j][2] =
              (constraints[i].layer-layer/2)*BAND_DIAM + layer_offset;
          }
        }
      }
    }
  }

  for (i = 0; i < n_constraints-1; i++) {
    if (constraints[i].radius) {
      vectorcp(constraints[i].crossings[constraints[i].n_crossings++],
               constraints[i].end_line);
    }
  }

  for (i = 0; i < n_constraints-1; i++) {
    if (constraints[i].radius) {
      calc_angles(band_type,&constraints[i],output);
    }
  }

  /*****************************************************
   * rotate everything back to whence it came
   * and put all parts back to original absolute
   * coordinates.
   *****************************************************/

  if ( ! ldraw_part) {
    fprintf(output,"0 SYNTH SYNTHESIZED BEGIN\n");
  }

  group_size = 0;

  /* now draw out the rubber band in terms of lines and arcs */
  for (i = 0; i < n_constraints; ) {
    if (constraints[i].radius != 0) {
      int j;
      for (j = i+1; j < n_constraints; j++) {

        if (constraints[j].radius) {
          draw_arc_line(
            band_type,
            &constraints[i],
            color,
            n_constraints > 1,
            output,
            ghost,
            group,
            &absolute,
            &constraints[first]);
          i = j;
          break;
        }
      }
      if (j >= n_constraints) {
        i++;
      }
    } else {
      i++;
    }
  }
  if (group) {
    fprintf(output,"0 GROUP %d %s\n",group_size,group);
  }
  if ( ! ldraw_part) {
    fprintf(output,"0 SYNTH SYNTHESIZED END\n");
  }
  return 0;

}

