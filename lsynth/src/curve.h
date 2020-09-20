/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef CURVE_H
#define CURVE_H

#ifdef _cplusplus
extern "C" {
#endif

int
synth_curve(
  part_t   *start,
  part_t   *end,
  part_t   *segments,
  int       n_segments,
  PRECISION attrib,
  FILE     *output);

PRECISION
hose_length(
  int       n_segments,
  part_t   *segments);
  
void
orient(
  part_t   *start,
  part_t   *end,
  int       n_segments,
  part_t   *segments);

void
orientq(
  part_t   *start,
  part_t   *end,
  int       n_segments,
  part_t   *segments);

#ifdef _cplusplus
};
#endif

#endif
