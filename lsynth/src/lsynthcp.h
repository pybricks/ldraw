/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef LSYNTH_H
#define LSYNTH_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mathlib.h"
#include "strings.h"

#define ACCY (1e-6)

typedef struct {
  char      type[128];
  PRECISION orient[3][3];
  PRECISION offset[3];
  PRECISION twist;
  int       attrib;
} part_t;

extern PRECISION max_bend;
extern PRECISION max_twist;
extern PRECISION band_res;
extern int group_size;
extern int ldraw_part;

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
  char            *type);

void list_products( void );

/************************************************************************
 *
 * Structures used to define the types of bands and hoses we can synthesize
 * and the constraints we can use to describe them.
 *
 **********************************************************************/

#if 1
// Somehow this broke RCX (STRETCH) cables in 3.1 beta g.
// Maybe what happened is not everything got recompiled...
#define STRETCH -1
#define FIXED   0
#define FIXED3  -2
#else
#define STRETCH 0
#define FIXED   1
#define FIXED3  2
#endif

#endif
