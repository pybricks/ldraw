void orient(
  part_t       *start,
  part_t       *end,
  int           n_segments,
  part_t       *segments)
{
  PRECISION start_up[3],end_up[3],up[3];
  PRECISION total_length = hose_length(n_segments,segments);
  PRECISION cur_length;
  int i;

  start_up[0] = 1; // Why is X the up vector?
  start_up[1] = 0; // The base constraint points up -Y,
  start_up[2] = 0; // And the up fin lies on -Z.
  rotate_point(start_up,start->orient);

  end_up[0] = 1;
  end_up[1] = 0;
  end_up[2] = 0;
  rotate_point(end_up,end->orient);

  /* Up vector controls the twist
   *
   * Interpolate the up vector based on start up vector, and
   * end up vector, and how far we are down the hose's total length
   */

  for (i = 0; i < n_segments-1; i++) {
    PRECISION r;
    PRECISION front[3];
    PRECISION t[3];

    cur_length = hose_length(i,segments);

    cur_length /= total_length;

    //Interpolate up vector from start to end based on length traversed.
    up[0] = start_up[0]*(1-cur_length) + end_up[0]*cur_length;
    up[1] = start_up[1]*(1-cur_length) + end_up[1]*cur_length;
    up[2] = start_up[2]*(1-cur_length) + end_up[2]*cur_length;

    r = sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
    up[0] /= r;
    up[1] /= r;
    up[2] /= r;

    // Get the current forward direction vector.
    front[0] = segments[i+1].offset[0] - segments[i].offset[0];
    front[1] = segments[i+1].offset[1] - segments[i].offset[1];
    front[2] = segments[i+1].offset[2] - segments[i].offset[2];

    r = sqrt(front[0]*front[0] + front[1]*front[1] + front[2]*front[2]);

    front[0] /= r;
    front[1] /= r;
    front[2] /= r;

    //Use cross products to create 2 Euler rotation axes?  gimble lock?
    mult_point(t,front,up); // side
    mult_point(up,t,front); // side * front

    //Take cross product (Up x Front) and make the rotation matrix from it.
    make_rotation_pp(segments[i].orient,up,front);
  }
}

    // ----------------------------------------------------------------
    // We have a problem.  

    // Interpolated up and forward direction vector may not be perpendicular.

    // Using quaternions
    // ----------------
    // Using forward vectors S (start) and V (front):
    // Calculate turn axis from cross product.  (S x V) / |(S x V)|
    // Then get the total turn angle from dot product.  arccos((S . V)/ |S|*|V|)
    // Make quaternion Q from axis and angle.  
    // Make rotation R matrix from quaternion.
    // Multiply start matrix by new matrix to get orient for segment.

    // Now for the spin:
    // ----------------
    // Call whatever spin is built into S a start spin of 0.
    // On the last segment above V = E.  Keep the last rotation matrix R.
    // Using up vectors s (start) and e (end):
    // Rotate s by R to get u.
    // Use dot product of u and e to get delta spin angle from s to e.
    // Interpolate the delta spin along the length of the hose.
    // We can apply the interpolated spin as a rotation around Y axis before 
    // applying the orient matrix of each segment.  
    // (Just add it to the per segment spin we already have).
    // ----------------------------------------------------------------
    

