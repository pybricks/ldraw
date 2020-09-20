/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague
 */

#include "lsynthcp.h"
#include "hose.h"
#include "curve.h"
#include "mathlib.h"

/*
 * a 1x1 brick is 20 LDU wide and 24 LDU high
 *
 * hose length = 14 brick widths long = 280 LDU
 * number of ribs = 45
 * 6.2 LDU per rib
 *
 */

void orient(
  int     n_segments,
  part_t *segments)
{
  int i;

  for (i = 0; i < n_segments-1; i++) {

      PRECISION dx, dy, dz;
      PRECISION orient[3][3];
      PRECISION L1,L2;

      // determine the orientation of the part in the XY plane
      // WARNING!  This is bad for twist if the hose mostly bends in XZ or YZ plane.
      // Start with the orient of the 1st constraint (the fin), and adjust from there?
      // Quaternion roll interpolation?  
      // Convert constraint orientations to quaternions and slerp the points between?

      dx = segments[i+1].offset[0] - segments[i].offset[0];
      dy = segments[i+1].offset[1] - segments[i].offset[1];
      dz = segments[i+1].offset[2] - segments[i].offset[2];

      L1 = sqrt(dx*dx + dy*dy);
      L2 = sqrt(dx*dx + dy*dy + dz*dz);
      if ((L1 == 0) || (L2 == 0)){ // Avoid divide by 0.
        segments[i].orient[0][0] = 1;
        segments[i].orient[1][0] = 0;
        segments[i].orient[2][0] = 0;
        segments[i].orient[0][1] = 0;
        segments[i].orient[1][1] = 0;
        segments[i].orient[2][1] =-1;
        segments[i].orient[0][2] = 0;
        segments[i].orient[1][2] = 1;
        segments[i].orient[2][2] = 0;
      } else {
        segments[i].orient[0][0] =  dy/L1;  //  cos
        segments[i].orient[1][0] = -dx/L1;  //  sin
        segments[i].orient[2][0] =  0;
        segments[i].orient[0][1] =  dx/L2;  // -sin
        segments[i].orient[1][1] =  dy/L2;  //  cos
        segments[i].orient[2][1] =  dz/L2;
        segments[i].orient[0][2] = -dx*dz/(L1*L2);
        segments[i].orient[1][2] = -dy*dz/(L1*L2);
        segments[i].orient[2][2] =  L1/L2;
      }
#if 1
      if ((i== 0) && (n_segments < 100)) 
	printf("Orient seg[%d](%.1f, %.1f, %.1f) = \n  (%.1f, %.1f, %.1f) (%.1f, %.1f, %.1f) (%.1f, %.1f, %.1f)\n", 
	       i, dx,dy,dz,
	       segments[i].orient[0][0],
	       segments[i].orient[1][0],
	       segments[i].orient[2][0],
	       segments[i].orient[0][1],
	       segments[i].orient[1][1],
	       segments[i].orient[2][1],
	       segments[i].orient[0][2],
	       segments[i].orient[1][2],
	       segments[i].orient[2][2]);
#endif

      // Try the more robust? (but ugly) reorient code from ldglite hoser fn.
      // This does not seem to do any better. (also operates in XY plane?)
      if (n_segments < 100)
      {
	double Xaxisrotation,Yaxisrotation,pi,pidiv2;
	pi=acos(-1);
	pidiv2=pi/2;
	
	if (dx<0){
	  if (dy<0){
	    Xaxisrotation=-(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	    if (dz<0) {
	      Yaxisrotation=(atan((dx)/(dz)));}
	    else if (dz>0) {
	      Yaxisrotation=(-pi+(atan((dx)/(dz))));}
	    else {
	      Yaxisrotation=pidiv2;}
	  }
	  else if (dy>0){
	    if (dz<0) {
	      Xaxisrotation=pi-(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	      Yaxisrotation=(atan((dx)/(dz)));}
	    else if (dz>0) {
	      Xaxisrotation=pi-(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	      Yaxisrotation=(-pi+(atan((dx)/(dz))));}
	    else {
	      Xaxisrotation=(-pi+(atan((dx)/(dy))));
	      Yaxisrotation=pidiv2;}
	  }
	  else {// if (dy==0)
	    Xaxisrotation=pidiv2;
	    if (dz<0) {
	      Yaxisrotation=(pi+(atan((dx)/(dz))));}
	    else if (dz>0) {
	      Yaxisrotation=(pi+(atan((dx)/(dz))));}
	    else {
	      Yaxisrotation=pidiv2;}
	  }
	}
	//***********************************************************
	else if (dx>0){
	  if (dy<0){
	  Xaxisrotation=-(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	  if (dz<0) {
	    Yaxisrotation=(atan((dx)/(dz)));}
	  else if (dz>0) {
	    Yaxisrotation=(-pi+(atan((dx)/(dz))));}
	  else {
	    Yaxisrotation=pidiv2;}
	}
	else if (dy>0){
	  Xaxisrotation=pi-(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	  if (dz<0) {
	    Yaxisrotation=(atan((dx)/(dz)));}
	  else if (dz>0) {
	    Yaxisrotation=(-pi+(atan((dx)/(dz))));}
	  else {
	    Yaxisrotation=pidiv2;}
	}
	else { // if (dy==0)
	  Xaxisrotation=-pidiv2;
	  if (dz<0) {
	    Yaxisrotation=(atan((dx)/(dz)));}
	  else if (dz>0) {
	    Yaxisrotation=(atan((dx)/(dz)));}
	  else {
	    Yaxisrotation=-pidiv2;}
	}
      }
      //***********************************************************
      else{ // if (dx==0)
	if (dy<0){
	  if (dz<0) {
	    Xaxisrotation=-(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	    Yaxisrotation=0;}
	  else if (dz>0) {
	    Xaxisrotation=(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	    Yaxisrotation=0;}
	  else {
	    Xaxisrotation=0;
	    Yaxisrotation=pi;}
	}
	else if (dy>0){
	  if (dz<0) {
	    Xaxisrotation=pi-(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	    Yaxisrotation=0;}
	  else if (dz>0) {
	    Xaxisrotation=-pi+(atan(sqrt(((dx)*(dx))+((dz)*(dz)))/(dy)));
	    Yaxisrotation=0;}
	  else {
	    Xaxisrotation=pi;
	    Yaxisrotation=pi;}
	}
	else { // if (dy==0)
	  if (dz<0) {
	    Xaxisrotation=pidiv2;
	    Yaxisrotation=0;}
	  else if (dz>0) {
	    Xaxisrotation=-pidiv2;
	    Yaxisrotation=0;}
	  else {
	    Xaxisrotation=pi;
	    Yaxisrotation=pi;}
	}
      }
	
	// Turn the parts 180 degrees to make them face the same way.
	if (Xaxisrotation < 0)
	  Xaxisrotation += pi; 
	else
	  Xaxisrotation -= pi; 

	segments[i].orient[0][0] = cos(Yaxisrotation);
	segments[i].orient[0][1] = (sin(Xaxisrotation)*sin(Yaxisrotation));
	segments[i].orient[0][2] = (cos(Xaxisrotation)*sin(Yaxisrotation));
	segments[i].orient[1][0] = 0;
	segments[i].orient[1][1] = cos(Xaxisrotation);
	segments[i].orient[1][2] = -sin(Xaxisrotation);
	segments[i].orient[2][0] = -sin(Yaxisrotation);
	segments[i].orient[2][1] = (sin(Xaxisrotation)*cos(Yaxisrotation));
	segments[i].orient[2][2] = (cos(Xaxisrotation)*cos(Yaxisrotation));

	if ((i== 0) && (n_segments < 100)) 
	  printf("ORIENT seg[%d](%.1f, %.1f, %.1f) = \n  (%.1f, %.1f, %.1f) (%.1f, %.1f, %.1f) (%.1f, %.1f, %.1f)\n", 
	       i, dx,dy,dz,
	       segments[i].orient[0][0],
	       segments[i].orient[1][0],
	       segments[i].orient[2][0],
	       segments[i].orient[0][1],
	       segments[i].orient[1][1],
	       segments[i].orient[2][1],
	       segments[i].orient[0][2],
	       segments[i].orient[1][2],
	       segments[i].orient[2][2]);
      }
  }
}


#if 0
      // quaternion(a,b,c,s) is just a direction vector with a spin(roll) angle.
      // ie. quaternion(x,y,z,spin)

      // convert rotation matrix m[16] to quaternion(a,b,c,s).
      // Calculate the trace of the matrix T from the equation:
      T = 1 + m[0] + m[5] + m[10];
      if ( T > 0.00000001 )
      {
        S = sqrt(T) * 2;
        a = ( m[9] - m[6] ) / S;
        b = ( m[2] - m[8] ) / S;
        c = ( m[4] - m[1] ) / S;
        s = 0.25 * S;
      }
      // If the trace of the matrix is equal to zero then identify
      // which major diagonal element has the greatest value.
      // Depending on this, calculate the following:
      else if ( m[0] > m[5] && m[0] > m[10] )  {	// Column 0: 
        S  = sqrt( 1.0 + m[0] - m[5] - m[10] ) * 2;
        a = 0.25 * S;
        b = (m[4] + m[1] ) / S;
        c = (m[2] + m[8] ) / S;
        s = (m[9] - m[6] ) / S;
      } else if ( m[5] > m[10] ) {			// Column 1: 
        S  = sqrt( 1.0 + m[5] - m[0] - m[10] ) * 2;
        a = (m[4] + m[1] ) / S;
        b = 0.25 * S;
        c = (m[9] + m[6] ) / S;
        s = (m[2] - m[8] ) / S;
      } else {						// Column 2:
        S  = sqrt( 1.0 + m[10] - m[0] - m[5] ) * 2;
        a = (m[2] + m[8] ) / S;
        b = (m[9] + m[6] ) / S;
        c = 0.25 * S;
        s = (m[4] - m[1] ) / S;
      }

      // convert quaternion(a,b,c,s) to a rotation matrix m[16].
      m[0]  = 1 - 2*b*b-2*c*c;
      m[1]  = 2*a*b - 2*s*c;
      m[2]  = 2*a*c + 2*s*b;
      m[4]  = 2*a*b+2*s*c;
      m[5]  = 1 - 2*a*a - 2*c*c;
      m[6]  = 2*b*c - 2*s*a;
      m[8]  = 2*a*c - 2*s*b;
      m[9]  = 2*b*c + 2*s*a;
      m[10] = 1 - 2*a*a - 2*b*b;
      m[3]  = m[7] = m[11] = m[12] = m[13] = m[14] = 0;
      m[15] = 1;
#endif



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
  PRECISION up[3];

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

  orient(n_segments, segments);

  return 0; /* it all worked */
}

