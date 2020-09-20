/*
 * This is the LDRAW parts synthesis library.
 * By Kevin Clague and Don Heyse
 */

#include <math.h>

#include "lsynthcp.h"
#include "hose.h"
#include "curve.h"
#include "mathlib.h"

#define PI 2*atan2(1,0)

/*
 * 0 SYNTH BEGIN DEFINE HOSE <fill> PNEUMATIC_HOSE "descr" <diameter> <stiffness> <twist>
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 0 SYNTH END
 */

int           n_hose_types = 0;
hose_attrib_t   hose_types[64];

#define N_HOSE_TYPES n_hose_types

/* In hoses, the attrib field in constraints, indicates that
 * LSynth should turn the final constraint around to get everything
 * to work correctly (e.g. flex-axle ends).
 */

/*
 * 0 SYNTH BEGIN DEFINE HOSE CONSTRAINTS
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 1 <length> a b c  d e f  g h i  j k l <part>
 * 0 SYNTH END
 */

int    n_hose_constraints = 0;
part_t   hose_constraints[128];

#define N_HOSE_CONSTRAINTS n_hose_constraints

void
list_hose_types(void)
{
  int i;

  printf("\n\nHose like synthesizable parts\n");
  for (i = 0; i < N_HOSE_TYPES; i++) {
    printf("  %-20s %s\n",hose_types[i].type, hose_types[i].descr);
  }
}

void
list_hose_constraints(void)
{
  int i;

  printf("\n\nHose constraints\n");
  for (i = 0; i < N_HOSE_CONSTRAINTS; i++) {
    printf("    %11s\n",hose_constraints[i].type);
  }
}

void
hose_ini(void)
{
  int i;

  for (i = 0; i < N_HOSE_TYPES; i++) {
    printf("%-20s = SYNTH BEGIN %s 16\n",hose_types[i].type, hose_types[i].type);
  }
}

int
ishosetype(char *type)
{
  int i;

  for (i = 0; i < N_HOSE_TYPES; i++) {
    if (strncasecmp(hose_types[i].type,type,strlen(hose_types[i].type)) == 0) {
      return 1;
    }
  }
  return 0;
}
// casecmp
int
ishoseconstraint(char *type)
{
  int i;

  for (i = 0; i < N_HOSE_CONSTRAINTS;i++) {
    if (strcasecmp(hose_constraints[i].type,type) == 0) {
      return 1;
    }
  }
  return 0;
}

PRECISION
line_angle(
  PRECISION va[3],
  PRECISION vb[3])
{
  PRECISION denom;
  PRECISION theta;

  denom = vectorlen(va)*vectorlen(vb);

  if (fabs(denom) < 1e-9) {
    //printf("line angle calculation gives us length of zero\n");
    return 0;
  }

  theta = (va[0]*vb[0]+va[1]*vb[1]+va[2]*vb[2])/denom;
  if (theta >= 1 || theta <= -1) {
    theta = 0;
  } else {
    theta = acos(theta);
  }
  return theta;
}

PRECISION
line_angle3(
  int a,
  int b,
  part_t *segments)
{
  PRECISION va[3]; /* line A */
  PRECISION vb[3]; /* line B */
  PRECISION denom;
  PRECISION theta;

  /* we get the line A as point[a] to point[a+1] */
  /* and line B as point[b] to point[b+1] */

  vectorsub3(va,segments[a+1].offset,segments[a].offset);
  vectorsub3(vb,segments[b+1].offset,segments[b].offset);

  return line_angle(va,vb);
}

/*
 * merge adjacent points until we either experience too much
 * bend, or too much twist
 */

int
merge_segments_angular(
  part_t    *start,
  part_t    *end,
  part_t    *segments,
  int       *n_segments,
  PRECISION  max_bend,
  PRECISION  max_twist,
  FILE      *output)
{
  int a,b;
  int n;
  PRECISION theta1,theta2;
  PRECISION total_length, cur_length;
  PRECISION start_up[3] = { 1, 0, 0 };
  PRECISION end_up[3]   = { 1, 0, 0 };
  PRECISION cur_up[3]   = { 1, 0, 0 };
  PRECISION next_up[3];
  PRECISION len[3], sub_len;

  vectorrot(start_up,start->orient);
  vectorrot(cur_up,  start->orient);
  vectorrot(end_up,  end->orient);
  total_length = hose_length(*n_segments, segments);
  cur_length = 0;

  a = 0; b = 1;
  n = 0;

  do {
    vectorsub3(len,segments[a].offset,segments[b].offset);
    sub_len = vectorlen(len);
    if (sub_len < 1e-9) {
      b++;
    }
  } while (sub_len < 1e-9);
  
  while (b < *n_segments) {
    PRECISION len[3],normalized;

    if (b < *n_segments - 1) {
      vectorsub3(len,segments[b+1].offset,segments[b].offset);
      cur_length += vectorlen(len);
    } else {
      cur_length += (total_length - cur_length)/2;
    }

    if (cur_length != 0) {
      PRECISION left;

      normalized = cur_length/total_length;
      left = 1 - normalized;

      next_up[0] = start_up[0]*left + end_up[0]*normalized;
      next_up[1] = start_up[1]*left + end_up[1]*normalized;
      next_up[2] = start_up[2]*left + end_up[2]*normalized;

      normalized = next_up[0]*next_up[0] +
                   next_up[1]*next_up[1] +
                   next_up[2]*next_up[2];
      normalized = sqrt(normalized);
      next_up[0] /= normalized;
      next_up[1] /= normalized;
      next_up[2] /= normalized;

      theta2 = line_angle(cur_up,next_up);
      theta1 = line_angle3(a,b,segments);

      if (theta1 < max_bend && theta2 < max_twist) {
        b++;
      } else {
        segments[n++] = segments[a++];
        a = b++;
        cur_up[0] = next_up[0];
        cur_up[1] = next_up[1];
        cur_up[2] = next_up[2];
      }
    } else {
      b++;
    }
  }
  if (n <= 2) {
    segments[1] = segments[*n_segments-1];
    n = 2;
  } else if (b - a > 1) {
    //n--; // -= 2;
    segments[n-1] = segments[a];
  } else {
    n--;
  }
  *n_segments = n;
  return 0;
}

int
merge_segments_length(
  part_t    *segments,
  int       *n_segments,
  PRECISION  max,
  FILE      *output)
{
  int a,b;
  int n;
  PRECISION d[3],l;

  a = 0; b = 1;
  n = 1;

  while (b < *n_segments) {

    vectorsub3(d,segments[a].offset,segments[b].offset);
    l = vectorlen(d);

    if (l + 0.5 < max) {
      b++;
    } else {
      a = b;
      segments[n++] = segments[b++];
    }
  }
  if (n < 2) {
    segments[1] = segments[*n_segments-1];
    n = 2;
  } else if (b - a > 1) {
    //n--; // -= 2;
    segments[n-1] = segments[a];
  } else {
    n--;
  }
  *n_segments = n;
  return 0;
}

int
merge_segments_count(
  hose_attrib_t  *hose,
  part_t    *start,
  part_t    *end,
  part_t    *segments,
  int       *n_segments,
  int       count,
  FILE      *output)
{
  int n, i;
  PRECISION d[3],l;
  PRECISION len, lenS, lenM, lenE;

  // Get the total length of the curve and divide by the expected segment count.
  len = 0;
  for (i = 0; i < *n_segments-1; i++) {
    vectorsub3(d,segments[i].offset,segments[i+1].offset);
    len += vectorlen(d);
  }
  printf("Total segment len = %.3f\n", len);

  // If S or E do not match the N parts subtract their lengths from total.
  if ((strcasecmp(hose->start.type, hose->mid.type) != 0) ||
      (strcasecmp(hose->end.type, hose->mid.type) != 0) ||
      (hose->start.attrib != hose->mid.attrib) ||
      (hose->end.attrib != hose->mid.attrib))
  {
    lenS = hose->start.attrib;
    lenE = hose->end.attrib;
    lenM = hose->mid.attrib;
    len = len - (lenS + lenE);
    printf("Net segment len = %.3f (S=%d, M=%d, E= %d)\n", len, lenS, lenM, lenE);
    len = len / (PRECISION)(count-1); //len /= (count);
  }
  else
  {
    len = len / (PRECISION)(count-1); //len /= (count);
    lenS = lenE = len;
  }
  printf("Merging %d segments to %d segments of len %.3f\n", *n_segments, count, len);

  // Break up the curve into count intervals of length len.
  l = 0;
  n = 1; // Keep the first point.
  // Find intermediate points.
  for (i = 0; i < *n_segments-1; i++) {
    vectorsub3(d,segments[i].offset,segments[i+1].offset);
    l += vectorlen(d);
    
    if ((l + 0.05) > (((n-1) * len) + lenS))
      segments[n++] = segments[i+1];

    //if (n >= count) break;
  }

  if (lenE != len) // If E did not match above, place it at the constraint.
  {
    segments[n-1] = segments[*n_segments-1]; // Use the last point twice?
  }
  segments[n++] = segments[*n_segments-1]; // Keep the last point.
  *n_segments = n;

  // NOTE: What I really need to do here is place the last point
  // segments[n-1].offset at the location of the end constraint + a lenE 
  // offset in the direction of the orientation of the end constraint.
  // otherwise the final vector is tiny and rounding can point it backwards.
  // We can see a backwards vector by checking for a negative dot product.

  // So, move the last point lenE from the start of the End constraint.
  // This should place it at the far point of the end constraint.
  d[0] = 0; d[1] = lenE; d[2] = 0;   // Create an offset vector of lenE
  d[1] *= -1;                        // along the -Y axis.   Reorient it
  vectorrot(d,end->orient);          // along the end constraint axis, and
  vectoradd(segments[n-1].offset,d); // add it to the end constraint origin.

  if (0) // Debug printouts
  {
    extern PRECISION dotprod(PRECISION a[3], PRECISION b[3]); // from curve.c

    PRECISION last, next;
    PRECISION dn[3],dp;
    PRECISION m2[3][3];

    matrixcp(m2,end->orient); 
    printf("  E =(%g, %g, %g, %g, %g, %g, %g, %g, %g)\n", m2[0][0], m2[0][1], m2[0][2],
	   m2[1][0], m2[1][1], m2[1][2], m2[2][0], m2[2][1], m2[2][2]);
    printf("  d=(%g, %g, %g)", d[0], d[1], d[2]);

    vectorsub3(d,segments[n-1].offset,segments[n-2].offset);
    last = vectorlen(d);
    vectorsub3(dn,segments[n-2].offset,segments[n-3].offset);
    next = vectorlen(dn);
    dp = dotprod(dn,d);
    printf("  last=%g, next=%g, dp=%g\n", last, next, dp);
  }

  printf("Produced %d points (%d segments)\n", *n_segments, *n_segments-1);

  // Reorient the segments.  
  // Warning!  Can interact badly with twist if hose makes a dx/dz (dy=0) turn.
  // Also, I think this ignores the orientation of the start and end constraints.
  // The fin on the constraints should guide the orientation (or twist?) somehow.
  // orient(n,segments);
  return 0;
}

#define MAX_SEGMENTS 1024*8

part_t segments[MAX_SEGMENTS];

// Create a second list to combine all patches between constraints.
part_t seglist[MAX_SEGMENTS];

void
output_line(
  FILE           *output,
  int             ghost,
  char           *group,
  int             color,
  PRECISION       a,
  PRECISION       b,
  PRECISION       c,
  PRECISION       d,
  PRECISION       e,
  PRECISION       f,
  PRECISION       g,
  PRECISION       h,
  PRECISION       i,
  PRECISION       j,
  PRECISION       k,
  PRECISION       l,
  char            *type)
{
  if (group) {
    fprintf(output,"0 MLCAD BTG %s\n",group);
  }
  fprintf(output,"%s1 %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %s\n",
      ghost ? "0 GHOST " : "",
      color,
      a,b,c,d,e,f,g,h,i,j,k,l,
      type);
  group_size++;
  fflush(stdout);
}

/*
 * Twist
 *    cos(t) 0 sin(t)
 *         0 1 0
 *   -sin(t) 0 cos(t)
 *
 * We need to add orientation of the segments based on the orientation
 * of the constraint.  We must bring the constraint into the normalized
 * orientation of the constraint (pointing up at Y), and determine the
 * angle of the tab relative to the Y axis.
 */

void
render_hose_segment(
  hose_attrib_t  *hose,
  int             ghost,
  char           *group,
  int            *group_size,
  int             color,
  part_t         *segments,
  int             n_segments,
  PRECISION      *total_twist,
  int             first,
  int             last,
  FILE           *output,
  part_t         *constraint)
{
  int i,j,k;
  PRECISION pi = 2*atan2(1,0);
  PRECISION m1[3][3];
  PRECISION m2[3][3];
  PRECISION offset[3];
  char     *type;
  int       gs = *group_size;

  for (i = 0; i <  n_segments-1; i++) {
    PRECISION tx,ty,tz,l,theta;

    if (hose->fill != STRETCH) {
      l = 1;
    } else {
      PRECISION d[3];

      vectorsub3(d,segments[i+1].offset,segments[i].offset);

      // the length of the segment is the distance between the
      // two points

      l = vectorlen(d);

      // plus the length needed to make sure that the outer edges
      // of the cross sections touch

      if (n_segments > 1) {
        PRECISION phi;

        phi = line_angle3(i+1,i,segments);
        phi = sin(phi)*(hose->diameter);

        l += phi + phi;
      }
    }

    m1[0][0] = 1;
    m1[0][1] = 0;
    m1[0][2] = 0;
    m1[1][0] = 0;
    m1[1][1] = l; // Stretch it by length of segment.
    m1[1][2] = 0;
    m1[2][0] = 0;
    m1[2][1] = 0;
    m1[2][2] = 1;

    // Default offset is nothing.
    offset[0] = 0; offset[1] = 0; offset[2] = 0;

    if (i == 0 && first) {
      type = hose->start.type;
      vectorcp(offset,hose->start.offset);
      if (hose->start.attrib != 0) // FIXED size start, do not stretch it. 
	matrixcp(m2,hose->start.orient);
      else  // Stretch it.
	matrixmult3(m2,hose->start.orient,m1);
    } else if (i == n_segments-2 && last) {
      type = hose->end.type;
      vectorcp(offset,hose->end.offset);
      if (hose->end.attrib != 0) // FIXED size end, do not stretch it. 
	matrixcp(m2,hose->end.orient);
      else // Stretch it.
	matrixmult3(m2,hose->end.orient,m1);
    } else if ((i & 0x01) && (strlen(hose->alt.type) != 0)) {
      if (i == 1) printf("ALT = %s\n", hose->alt.type);
      type = hose->alt.type;
      vectorcp(offset,hose->alt.offset);
      matrixmult3(m2,hose->alt.orient,m1);
    } else {
      type = hose->mid.type;
      vectorcp(offset,hose->mid.offset);
      matrixmult3(m2,hose->mid.orient,m1);
    }

    /*
     * twist helps with string or chains
     */

    m1[0][0] = 1;
    m1[0][1] = 0;
    m1[0][2] = 0;
    m1[1][0] = 0;
    m1[1][1] = 1;
    m1[1][2] = 0;
    m1[2][0] = 0;
    m1[2][1] = 0;
    m1[2][2] = 1;

    if (hose->fill != STRETCH) {
      PRECISION angle;
      *total_twist += hose->twist;       // One twist before calculating angle.
#define ORIENT_FN_FIXED_FOR_XZ_AND_YZ_CURVES
#ifdef ORIENT_FN_FIXED_FOR_XZ_AND_YZ_CURVES
      // For N FIXED segments start with a full twist (not a half twist).
      if (hose->fill > FIXED) {
	angle = *total_twist * pi / 180; // Calculate angle after the twist.
      }else 
#endif
      {
	angle = *total_twist * pi / 360; // Not (pi/180) because we double the twist.
	*total_twist += hose->twist;     // Another twist after calculating angle.
	// NOTE: 
	//       Breaking the twist up this way draws the Minifig chain with
	//       the links twisted 45 degrees from the start post (instead of 90).
	//       This looks OK because the real chains sometimes rest that way.
	//       It also avoids the problem with the orient() fn.
	//       Chains which turn in the XZ or YZ planes may rotate 45 degrees
	//       the other way because orient() only works in the XY plane.
	//       But that's better than the 90 degrees wrong without this hack.
      }
      m1[0][0] =   cos(angle);
      m1[0][2] =   sin(angle);
      m1[2][0] =  -sin(angle);
      m1[2][2] =   cos(angle);
    }
    matrixmult(m1,m2);
    matrixmult3(m2,segments[i].orient,m1);

    // NOTE: We handle hose->(start,mid,end).orient here, but we do nothing
    //       with the hose->(start,mid,end).offset.
    // I really think we need to use this to fix the minifig chain link, the 
    // origin of which is not centered.  Instead its ~4Y over toward one end.
#if 0
    offset[0] = 0; offset[1] = 0; offset[2] = 0;
    vectoradd(offset,hose->mid.offset); // Should also consider first and last.
#else
    // Get offset (first, mid, or last) from lsynth.mpd file above.
    // That way it matches the first, mid, or last orient from lsynth.mpd.
#endif
    vectorrot(offset,m2);
    vectoradd(segments[i].offset,offset);

    // FIXME: I would expect to have to flip the array along the diagonal
    // like we have to do in input.

    output_line(output,ghost,group,color,
      segments[i].offset[0], segments[i].offset[1], segments[i].offset[2],
      m2[0][0], m2[0][1], m2[0][2],
      m2[1][0], m2[1][1], m2[1][2],
      m2[2][0], m2[2][1], m2[2][2],
      type);
    if (group) {
      gs++;
    }
  }
  *group_size = gs;
}

void
adjust_constraint(
  part_t *part,
  part_t *orig,
  int     last)
{
  int i;
  PRECISION m[3][3];

  *part = *orig;

  for (i = 0; i < N_HOSE_CONSTRAINTS; i++) {
    if (strcasecmp(part->type,hose_constraints[i].type) == 0) {

      // adjust the constraints offset via hose_constraint
      vectorcp(part->offset,hose_constraints[i].offset);
      vectorrot(part->offset,hose_constraints[i].orient);
      vectoradd(part->offset,orig->offset);

      // compensate part orient via hose_constraint
      matrixcp(part->orient,hose_constraints[i].orient);
      matrixmult(part->orient,orig->orient);

      if (hose_constraints[i].attrib && last) {
        matrixcp(m,part->orient);
        matrixneg(part->orient,m);
      }
      break;
    }
  }
  return;
}

/*
 * a 1x1 brick is 20 LDU wide and 24 LDU high
 *
 * hose length = 14 brick widths long = 280 LDU
 * number of ribs = 45
 * 6.2 LDU per rib
 *
 */

void
render_hose(
  hose_attrib_t  *hose,
  int             n_constraints,
  part_t         *constraints,
  PRECISION       bend_res,
  PRECISION       twist_res,
  int             ghost,
  char           *group,
  int             group_size,
  int             color,
  FILE *output)
{
  int       c, n_segments;
  part_t    mid_constraint;
  PRECISION total_twist = 0;

  if ( ! ldraw_part) {
    fprintf(output,"0 SYNTH SYNTHESIZED BEGIN\n");
  }

  // First and Last parts for STRETCH hose could be FIXED length.
  // Add up to two new constraints to handle this.  They're not used afterwards.
  // Just don't go over 128 constraints.
  if (hose->fill == STRETCH) {
    PRECISION offset[3];
    PRECISION l;

#ifdef DEBUGGING_HOSES
    printf("STRETCH = (%d, %d, %d)\n", 
	   hose->start.attrib, hose->mid.attrib, hose->end.attrib);
#endif

    l = hose->start.attrib;
    if (l != 0.0) {
      n_constraints++;
      for (c = n_constraints-1; c > 0; c--) {
	memcpy(&constraints[c], &constraints[c-1], sizeof(part_t));
      }

      offset[0] = 0; offset[1] = -l; offset[2] = 0;
      vectorrot(offset,constraints[0].orient);
      vectoradd(constraints[1].offset, offset);
    }

    l = hose->end.attrib;
    if (l != 0.0){
      n_constraints++;
      c = n_constraints-1;
      memcpy(&constraints[c], &constraints[c-1], sizeof(part_t));

      offset[0] = 0; offset[1] = l; offset[2] = 0;
      vectorrot(offset,constraints[n_constraints-1].orient);
      vectoradd(constraints[n_constraints-2].offset, offset);
    }
  }

  mid_constraint = constraints[0];

  for (c = 0; c < n_constraints - 1; c++) {
    part_t first,second;

    // reorient imperfectly oriented or displaced constraint types

    adjust_constraint(&first,&mid_constraint,0);

    // reorient imperflectly oriented or displaced constraint types

    adjust_constraint(&second, &constraints[c+1],c == n_constraints-2);

    n_segments = MAX_SEGMENTS;

    // For N FIXED segments break up segments[MAX_SEGMENTS] into chunks.
    if (hose->fill > FIXED)
    {
      // Break up seglist into fixed size chunks based on number of constraints.
      // We could do better by considering the actual length of each chunk...
      n_segments = MAX_SEGMENTS / (n_constraints - 1);
      //if (c == 0) printf("FIXED%d, N_constraints = %d, chunksize = %d\n", hose->fill, n_constraints, n_segments);
    }

    // create an oversampled curve

    if (hose->fill == FIXED) // Save room for end constraint point.
      synth_curve(&first,&second,segments,n_segments-1,hose->stiffness,output);
    else if (hose->fill == STRETCH) 
      synth_curve(&first,&second,segments,n_segments-1,hose->stiffness,output);
    else if (hose->fill > FIXED)
      synth_curve(&first,&second,segments,n_segments-1,hose->stiffness,output);
    else // Old way.  Overwrite last point with end constraint point.  Not good.
      synth_curve(&first,&second,segments,n_segments,hose->stiffness,output);

    // reduce oversampled curve to fixed length chunks, or segments limit
    // by angular resolution

    if (hose->fill == STRETCH) {
      // Make sure final segment matches second constraint
      vectorcp(segments[n_segments-1].offset,second.offset);
      merge_segments_angular(
        &first,
        &second,
        segments,
        &n_segments,
        bend_res,
        twist_res,
        output);
      // Make sure final segment matches second constraint
      vectorcp(segments[n_segments-1].offset,second.offset);
      // move normalized result back into its original orientation and position
      mid_constraint = constraints[c+1];
#ifdef DEBUGGING_HOSES
      printf("orient(N_SEGMENTS = %d)\n", n_segments);
#endif
      orientq(&first,&second,n_segments,segments); // With quaternions!
    }
    else if (hose->fill == FIXED) {
      // Make sure final segment point matches second constraint
      vectorcp(segments[n_segments-1].offset,second.offset);
#ifdef ADJUST_FINAL_FIXED_HOSE_END
      // Hmmm, how do we make a hose comprised of fixed length parts
      // reach exactly to the end constraint?
      // This is OK for ribbed hoses 
      // but not for string or chain, so we probably shouldn't do it.
      if (c == n_constraints-2) 
      {
	int i = n_segments;
	memcpy(seglist, segments, n_segments*sizeof(part_t));
	// Set i to how many segments we need to get near to the end.
	merge_segments_length(seglist,&i,hose->mid.attrib,output);
	// Squish an extra part into the last segment to make it reach the end.
	merge_segments_count(hose,&first,&second,segments,&n_segments,i,output); // Yuck!
	// Or, stretch the last part of the hose a bit to make it reach the end.
	// merge_segments_count(hose, segments,&n_segments,i-1,output); // Yuckier!
      }
      else
#endif
	merge_segments_length(segments,&n_segments,hose->mid.attrib,output);
      // move normalized result back into its original orientation and position
      mid_constraint = constraints[c+1];
      vectorcp(mid_constraint.offset,segments[n_segments-1].offset);
      orient(&first,&second,n_segments,segments);
      //orientq(&first,&second,n_segments,segments);
    }
    else { // For N fixed size chunks just copy into one big list, merge later.
      // Make sure final segment matches second constraint
      vectorcp(segments[n_segments-1].offset,second.offset);
      // move normalized result back into its original orientation and position
      mid_constraint = constraints[c+1];
      vectorcp(mid_constraint.offset,segments[n_segments-2].offset);
      memcpy(seglist+(n_segments*c), segments, n_segments*sizeof(part_t));
    }

    // output the result (if not FIXED number of segments)
    if (hose->fill <= FIXED)
    render_hose_segment(
      hose,
      ghost,
      group,
      &group_size,
      color,
      segments,n_segments,
      &total_twist,
      c ==0,
      c == n_constraints-2,
      output,
      &constraints[0]);
  }

  // output the result (if FIXED number of segments)
  if (hose->fill > FIXED) {
    part_t first,second;

    // Reorient imperfectly oriented or displaced constraint types.
    // NOTE: I really need to study the new orient code and see why it 
    // does not seem to work until after merge_segment_count() below.
    // It needs to orient based on ALL constraints, not just first and last.
    adjust_constraint(&first,&constraints[0],0);
    adjust_constraint(&second, &constraints[n_constraints-1],1);

    n_segments *= c;
    // Make sure final segment matches second constraint
    vectorcp(segments[n_segments-1].offset,second.offset);
    merge_segments_count(hose,&first,&second,seglist,&n_segments,hose->fill,output);
    //printf("Merged segments to %d segments of len %d\n", n_segments, hose->mid.attrib);
    //orient(&first,&second,n_segments,seglist);
    orientq(&first,&second,n_segments,seglist); // With quaternions!
    //printf("oriented %d\n",n_segments);
   
    render_hose_segment(
      hose,
      ghost,
      group,
      &group_size,
      color,
      seglist,n_segments,
      &total_twist,
      1, // First AND
      1, // Last part of the hose. (seglist = one piece hose)
      output,
      &constraints[0]);
      //printf("Total twist = %.1f (%.1f * %.1f) n = %d\n", total_twist, hose->twist, total_twist/hose->twist, n_segments);
  }

  if (group) {
    fprintf(output,"0 GROUP %d %s\n",group_size,group);
  }

  if ( ! ldraw_part) {
    fprintf(output,"0 SYNTH SYNTHESIZED END\n");
  }
}

int
synth_hose(
  char   *type,
  int     n_constraints,
  part_t *constraints,
  int     ghost,
  char   *group,
  int     group_size,
  int     color,
  FILE   *output)
{
  int i;

  for (i = 0; i < N_HOSE_TYPES; i++) {
    if (strcasecmp(hose_types[i].type,type) == 0) {
      render_hose(
        &hose_types[i],
        n_constraints,constraints,
        max_bend,
        max_twist,
        ghost,
        group,
        group_size,
        color,
        output);
      return 0;
    }
  }
  return 1;
}

