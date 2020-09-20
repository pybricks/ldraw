/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynthcp.h"
#include "hose.h"
#include "curve.h"
#include "mathlib.h"

static void
mult_point(PRECISION r[3], PRECISION lhs[3], PRECISION rhs[3])
{
  PRECISION tt;

  r[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
  r[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
  r[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];

  tt = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

  r[0] /= tt;
  r[1] /= tt;
  r[2] /= tt;
}

static void
make_rotation_pp(
  PRECISION r[3][3],
  PRECISION up[3],
  PRECISION front[3])
{
  PRECISION side[3];

  mult_point(side,front,up);

  r[0][0] = up[0];
  r[1][0] = up[1];
  r[2][0] = up[2];
  r[0][1] = front[0];
  r[1][1] = front[1];
  r[2][1] = front[2];
  r[0][2] = side[0];
  r[1][2] = side[1];
  r[2][2] = side[2];
}

static void
rotate_point(
  PRECISION r[3],
  PRECISION m[3][3])
{
  PRECISION t[3];

  t[0] = r[0]*m[0][0] + r[1]*m[0][1] + r[2]*m[0][2];
  t[1] = r[0]*m[1][0] + r[1]*m[1][1] + r[2]*m[1][2];
  t[2] = r[0]*m[2][0] + r[1]*m[2][1] + r[2]*m[2][2];

  r[0] = t[0];
  r[1] = t[1];
  r[2] = t[2];
}

void
orient2(
  part_t       *start,
  part_t       *end,
  int           n_segments,
  part_t       *segments)
{
  PRECISION up[3];
  int i;

  up[0] = 1;
  up[1] = 0;
  up[2] = 0;

  rotate_point(up,start->orient);
  printf("%5.2f %5.2f %5.2f\n",up[0],up[1],up[2]);

  for (i = 0; i < n_segments-1; i++) {
    PRECISION r;
    PRECISION front[3];
    PRECISION t[3];

    front[0] = segments[i+1].offset[0] - segments[i].offset[0];
    front[1] = segments[i+1].offset[1] - segments[i].offset[1];
    front[2] = segments[i+1].offset[2] - segments[i].offset[2];

    r = sqrt(front[0]*front[0] + front[1]*front[1] + front[2]*front[2]);

    front[0] /= r;
    front[1] /= r;
    front[2] /= r;

    mult_point(t,front,up); // side
    mult_point(up,t,front); // side * front

    make_rotation_pp(segments[i].orient,up,front);
  }
}

PRECISION hose_length(
  int           n_segments,
  part_t       *segments)
{
  PRECISION length = 0;
  int i;

  for (i = 0; i < n_segments-1; i++) {
    PRECISION l[3];

    vectorsub3(l,segments[i+1].offset,segments[i].offset);

    length += vectorlen(l);
  }
  return length;
}

void
orient(
  part_t       *start,
  part_t       *end,
  int           n_segments,
  part_t       *segments)
{
  PRECISION start_up[3],end_up[3],up[3];
  PRECISION total_length = hose_length(n_segments,segments);
  PRECISION cur_length;
  int i;

  start_up[0] = 1;
  start_up[1] = 0;
  start_up[2] = 0;
  rotate_point(start_up,start->orient);

  end_up[0] = 1;
  end_up[1] = 0;
  end_up[2] = 0;
  rotate_point(end_up,end->orient);

  /* Up vector controls the twist
   *
   * Interpolate the up vector based on start up vector, and
   * end up vector, and how far we are down the hose's total
   * length
   */

  for (i = 0; i < n_segments-1; i++) {
    PRECISION r;
    PRECISION front[3];
    PRECISION t[3];

    cur_length = hose_length(i,segments);

    cur_length /= total_length;

    up[0] = start_up[0]*(1-cur_length) + end_up[0]*cur_length;
    up[1] = start_up[1]*(1-cur_length) + end_up[1]*cur_length;
    up[2] = start_up[2]*(1-cur_length) + end_up[2]*cur_length;

    r = sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
    up[0] /= r;
    up[1] /= r;
    up[2] /= r;

    front[0] = segments[i+1].offset[0] - segments[i].offset[0];
    front[1] = segments[i+1].offset[1] - segments[i].offset[1];
    front[2] = segments[i+1].offset[2] - segments[i].offset[2];

    r = sqrt(front[0]*front[0] + front[1]*front[1] + front[2]*front[2]);

    front[0] /= r;
    front[1] /= r;
    front[2] /= r;

    mult_point(t,front,up); // side
    mult_point(up,t,front); // side * front

    make_rotation_pp(segments[i].orient,up,front);
  }
}

int
synth_curve(
  part_t       *start,
  part_t       *end,
  part_t       *segments,
  int           n_segments,
  PRECISION     attrib,
  FILE         *output)
{
  PRECISION vector[3];
  PRECISION start_speed_v[3];
  PRECISION stop_speed_v[3];
  PRECISION time,time2,i_time,i_time_sum,i_time_sum2,n_time;
  PRECISION x,x2,y,y2,z,z2;
  PRECISION ptp,ratio;
  PRECISION ptp_sum;
  int i,j,n;

  vector[0] = 0;
  vector[1] = attrib;
  vector[2] = 0;

  for (j = 0; j < 3; j++)
  {
     x = 0.0;
     for (i = 0; i < 3; i++) {
        x += start->orient[j][i] * vector[i];
     }
     start_speed_v[j] = x;
  }

  vector[0] = 0;
  vector[1] = attrib;
  vector[2] = 0;

  for (j = 0; j < 3; j++)
  {
     x = 0.0;
     for (i = 0; i < 3; i++) {
        x -= end->orient[j][i] * vector[i];
     }
     stop_speed_v[j] = x;
  }

  // stop_speed_v[2] = -stop_speed_v[2];

  for (i = 0; i < n_segments; i++) {
    time  = (PRECISION) i/ (PRECISION) n_segments;

    segments[i].offset[0] =
      (1 - time) * (1 - time) * (1 - time) * start->offset[0] +
      (1 - time) * (1 - time) * 3 * time * (start->offset[0] - start_speed_v[0]) +
      (1 - time) * 3 * time * time * (end->offset[0] - stop_speed_v[0]) +
            time * time * time * end->offset[0];

    segments[i].offset[1] =
      (1 - time) * (1 - time) * (1 - time) * start->offset[1] +
      (1 - time) * (1 - time) * 3 * time * (start->offset[1] - start_speed_v[1]) +
      (1 - time) * 3 * time * time * (end->offset[1] - stop_speed_v[1]) +
            time * time * time * end->offset[1];
/*
=(1-$A8)^3*D$4 + (1-$A8)^2*3*$A8*(D$4-D$3) +  (1-$A8)*3*$A8^2*(D$5-D$6) + $A8^3*D$5
 */
    segments[i].offset[2] =
      (1 - time) * (1 - time) * (1 - time) * start->offset[2] +
      (1 - time) * (1 - time) * 3 * time * (start->offset[2] - start_speed_v[2]) +
      (1 - time) * 3 * time * time * (end->offset[2] - stop_speed_v[2]) +
            time * time * time * end->offset[2];
  }

  ptp_sum = 0;

  for (i = 0; i < n_segments - 1; i++) {
    x = segments[i + 1].offset[0] - segments[i].offset[0];
    y = segments[i + 1].offset[1] - segments[i].offset[1];
    z = segments[i + 1].offset[2] - segments[i].offset[2];

    ptp = sqrt(x*x + y*y + z*z);
    ptp_sum += ptp;
  }

  /* G8 */

  i_time_sum = 0;
  for (i = 0; i < n_segments - 1; i++) {
    time  = (PRECISION) (i+1)/ (PRECISION) n_segments;

    x = segments[i + 1].offset[0] - segments[i].offset[0];
    y = segments[i + 1].offset[1] - segments[i].offset[1];
    z = segments[i + 1].offset[2] - segments[i].offset[2];

    ptp = sqrt(x*x + y*y + z*z);

    ratio = ptp*n_segments/ptp_sum;
    if (ratio == 0) {
      ratio = 1e-20;
    }
    i_time = 1.0/(n_segments*ratio);
    i_time_sum += i_time;
  }

  /* H, I, J */
  n_time = 0;
  i_time_sum2 = 0;
  for (i = 0; i < n_segments - 1; i++) {
    PRECISION foo;

    x = segments[i + 1].offset[0] - segments[i].offset[0];
    y = segments[i + 1].offset[1] - segments[i].offset[1];
    z = segments[i + 1].offset[2] - segments[i].offset[2];

    ptp = sqrt(x * x + y * y + z * z);  /* E */
    ratio = ptp*n_segments/ptp_sum;     /* F */
    if (ratio == 0) {
      ratio = 1e-20;
    }
    i_time = 1.0/(n_segments*ratio);    /* G */
    i_time_sum2 += i_time;

    foo = 1.0/n_segments;
    foo /= ratio;
    foo /= i_time_sum;

    segments[i].offset[0] =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[0] +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[0] - start_speed_v[0]) +
      (1 - n_time) * 3 * n_time * n_time * (end->offset[0] - stop_speed_v[0]) +
       n_time * n_time * n_time * end->offset[0];

    segments[i].offset[1] =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[1] +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[1] - start_speed_v[1]) +
      (1 - n_time) * 3 * n_time * n_time * (end->offset[1] - stop_speed_v[1]) +
       n_time * n_time * n_time * end->offset[1];

    segments[i].offset[2] =
      (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[2] +
      (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[2] - start_speed_v[2]) +
      (1 - n_time) * 3 * n_time * n_time * (end->offset[2] - stop_speed_v[2]) +
       n_time * n_time * n_time * end->offset[2];

    n_time += foo;  /* H */
  }

  segments[i].offset[0] =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[0] +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[0] - start_speed_v[0]) +
    (1 - n_time) * 3 * n_time * n_time * (end->offset[0] - stop_speed_v[0]) +
     n_time * n_time * n_time * end->offset[0];

  segments[i].offset[1] =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[1] +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[1] - start_speed_v[1]) +
    (1 - n_time) * 3 * n_time * n_time * (end->offset[1] - stop_speed_v[1]) +
     n_time * n_time * n_time * end->offset[1];

  segments[i].offset[2] =
    (1 - n_time) * (1 - n_time) * (1 - n_time) * start->offset[2] +
    (1 - n_time) * (1 - n_time) * 3 * n_time * (start->offset[2] - start_speed_v[2]) +
    (1 - n_time) * 3 * n_time * n_time * (end->offset[2] - stop_speed_v[2]) +
     n_time * n_time * n_time * end->offset[2];

  // orient(n_segments, segments);

  return 0; /* it all worked */
}





/***************************************************************/
// Quaternion q[4] is mostly a direction vector with a spin(roll) angle.
// ie. quaternion(x,y,z,spin)
// 
// Gotta compare quat FAQ conversions to ldglite conversions to
// see what to do with the weird ldraw coordinate system.  I think I 
// may have to reverse -y (and maybe z) to get right? handed system
// before using quaternions.  Make sure I copied all conversions from
// same place to avoid mixed systems.

/***************************************************************/

/***************************************************************/
void normalizequat(PRECISION q[4])
{
  PRECISION L;

  L = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);

  q[0] /= L;
  q[1] /= L;
  q[2] /= L;
  q[3] /= L;
}

/***************************************************************/
void normalize(PRECISION v[3])
{
  PRECISION L;

  L = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  v[0] /= L;
  v[1] /= L;
  v[2] /= L;
}

/***************************************************************/
// Convert axis v[3] and angle a to quaternion q[4].
/***************************************************************/
void axisangle2quat(
  PRECISION v[3], 
  PRECISION a,
  PRECISION q[4])
{
  PRECISION sina,cosa;
  
  a = a / 2.0;
  sina = sin(a);
  cosa = cos(a);
  q[0] = v[0] * sina;
  q[1] = v[1] * sina;
  q[2] = v[2] * sina;
  q[3] = cosa;
}

/***************************************************************/
// Convert rotation matrix m[3][3] to quaternion q[4].
/***************************************************************/
void matrix2quat(
  PRECISION m[3][3],
  PRECISION q[4])
{
  PRECISION T;
  PRECISION S;

  // Calculate the trace of the matrix T from the equation:
  T = 1 + m[0][0] + m[1][1] + m[2][2];
  if ( T > 0.00000001 )
  {
    S = sqrt(T) * 2;
    q[0] = ( m[2][1] - m[1][2] ) / S;
    q[1] = ( m[0][2] - m[2][0] ) / S;
    q[2] = ( m[1][0] - m[0][1] ) / S;
    q[3] = 0.25 * S;
  }
  // If the trace of the matrix is equal to zero then identify
  // which major diagonal element has the greatest value.
  // Depending on this, calculate the following:
  else if ( m[0][0] > m[1][1] && m[0][0] > m[2][2] )  {	// Column 0: 
    S  = sqrt( 1.0 + m[0][0] - m[1][1] - m[2][2] ) * 2;
    q[0] = 0.25 * S;
    q[1] = ( m[1][0] + m[0][1] ) / S;
    q[2] = ( m[0][2] + m[2][0] ) / S;
    q[3] = ( m[2][1] - m[1][2] ) / S;
  } else if ( m[1][1] > m[2][2] ) {			// Column 1: 
    S  = sqrt( 1.0 + m[1][1] - m[0][0] - m[2][2] ) * 2;
    q[0] = ( m[1][0] + m[0][1] ) / S;
    q[1] = 0.25 * S;
    q[2] = ( m[2][1] + m[1][2] ) / S;
    q[3] = ( m[0][2] - m[2][0] ) / S;
  } else {						// Column 2:
    S  = sqrt( 1.0 + m[2][2] - m[0][0] - m[1][1] ) * 2;
    q[0] = ( m[0][2] + m[2][0] ) / S;
    q[1] = ( m[2][1] + m[1][2] ) / S;
    q[2] = 0.25 * S;
    q[3] = ( m[1][0] - m[0][1] ) / S;
  }
}

/***************************************************************/
// Convert quaternion q[4] to rotation matrix m[3][3].
/***************************************************************/
void quat2matrix(
  PRECISION q[4],
  PRECISION m[3][3])
{
  PRECISION a,b,c,s;

  normalizequat( q );

  a = q[0];
  b = q[1];
  c = q[2];
  s = q[3];

  m[0][0]  = 1 - 2*b*b-2*c*c;
  m[0][1]  = 2*a*b - 2*s*c;
  m[0][2]  = 2*a*c + 2*s*b;
  m[1][0]  = 2*a*b+2*s*c;
  m[1][1]  = 1 - 2*a*a - 2*c*c;
  m[1][2]  = 2*b*c - 2*s*a;
  m[2][0]  = 2*a*c - 2*s*b;
  m[2][1]  = 2*b*c + 2*s*a;
  m[2][2] = 1 - 2*a*a - 2*b*b;
}

/***************************************************************/
// Convert quaternion q[4] to a rotation axis v[3] and angle a.
/***************************************************************/
 
PRECISION quat2axisangle(
  PRECISION q[4],
  PRECISION v[3],
  PRECISION *a)
{
  PRECISION cos_angle, sin_angle;

  normalizequat( q );

  cos_angle  = q[3];
  *a          = acos( cos_angle ) * 2;
  sin_angle  = sqrt( 1.0 - cos_angle * cos_angle );

  if ( fabs( sin_angle ) < 0.0005 )
    sin_angle = 1;

  v[0] = q[0] / sin_angle;
  v[1] = q[1] / sin_angle;
  v[2] = q[2] / sin_angle;

  return *a;
}

/***************************************************************/
// Convert start, end constraints and n points to oriented segments
// using quaternion.
/***************************************************************/
void
orientq(
  part_t       *start,
  part_t       *end,
  int           n_segments,
  part_t       *segments)
{
  PRECISION start_up,end_up,up; // Up is now a twist angle, not a vector.
  PRECISION start_q[4],end_q[4],q[4];
  PRECISION v[3];
  PRECISION total_length = hose_length(n_segments,segments);
  PRECISION cur_length;
  int i;

  // Get the start and end twist angles (up vectors)
  matrix2quat(start->orient, start_q);
  matrix2quat(end->orient, end_q);

#define DEBUG_QUAT_MATH 1
#ifdef DEBUG_QUAT_MATH
 {
  PRECISION pi = 2*atan2(1,0);
  PRECISION degrees = 180.0 / pi;

  quat2axisangle(start_q, v, &start_up);
  quat2axisangle(end_q, v, &end_up);
  printf("Start_q = (%.2f, %.2f, %.2f, %.2f) up = %.2f degrees\n",
	 start_q[0],start_q[1],start_q[2],start_q[3], start_up*degrees);
 }
#else
  start_up = acos(start_q[3]) * 2.0;
  end_up = acos(end_q[3]) * 2.0;
#endif

  /* Up angle controls the twist
   *
   * Interpolate the up angle based on start up angle, end up angle, 
   * and how far we are down the hose's total length.
   */

  for (i = 1; i < n_segments; i++) {
    PRECISION v0[3];

    cur_length = hose_length(i,segments);
    cur_length /= total_length;
    up = start_up*(1-cur_length) + end_up*cur_length;

    // Get the axis before segment i.
    v0[0] = segments[i].offset[0] - segments[i-1].offset[0];
    v0[1] = segments[i].offset[1] - segments[i-1].offset[1];
    v0[2] = segments[i].offset[2] - segments[i-1].offset[2];
    normalize(v0);

    // Get the axis after segment i.
    v[0] = segments[i+1].offset[0] - segments[i].offset[0];
    v[1] = segments[i+1].offset[1] - segments[i].offset[1];
    v[2] = segments[i+1].offset[2] - segments[i].offset[2];
    normalize(v);

    // Tangent at segment i is the average of the before and after axis.
    v[0] += v0[0];
    v[1] += v0[1]; 
    v[2] += v0[2];
    normalize(v);
    
    // Convert to a quaternion, and convert that to a rotation matrix.
    axisangle2quat(v, up, q);
    normalizequat( q );
    quat2matrix(q, segments[i].orient);

#ifdef DEBUG_QUAT_MATH
    {
      PRECISION pi = 2*atan2(1,0);
      PRECISION degrees = 180.0 / pi;
      
      printf("  V[%d] = (%.2fx, %.2fy, %.2fz, %.2f up) => Q(%.2f, %.2f, %.2f, %.2f)\n",
	     i,v[0],v[1],v[2], up*degrees, q[0],q[1],q[2],q[3]);
    }
#endif
  }
#ifdef DEBUG_QUAT_MATH
  {
    PRECISION pi = 2*atan2(1,0);
    PRECISION degrees = 180.0 / pi;
      
    printf("End_q = (%.2f, %.2f, %.2f, %.2f) up = %.2f degrees\n",
	   end_q[0],end_q[1],end_q[1],end_q[1], end_up*degrees);
  }
#endif
}

