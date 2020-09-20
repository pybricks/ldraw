/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef BAND_H
#define BAND_H

#ifdef _cplusplus
extern "C" {
#endif

typedef struct {
  part_t    part;
  PRECISION radius;
  int       inside;
  int       was_cross;
  int       cross;
  PRECISION start_line[3];
  PRECISION end_line[3];
  PRECISION crossings[8][3];
  int       n_crossings;
  int       layer;

  PRECISION start_angle[3];
  PRECISION end_angle[3];
  int       n_steps;
  PRECISION s_angle;

  int       band_constraint_n;
} LSL_band_constraint;

typedef struct {
  char      type[256];   // name of thing being specified (e.g. RUBBER_BAND)
  char      descr[256];
  int       fill;        // method for synthesizing
                         // FIXED - chain and tread are composed of fixed
                         //         length parts
                         // FIXED3 - special case for rubber tread.
                         // STRETCH - for rubber bands
  PRECISION scale;       // convert LDUs to number of parts
  PRECISION thresh;      // used to compute number of steps on arc
  int       pulley;      // initially for STRING on pulleys
  part_t    tangent;     // the part type used for tangent
  part_t    arc;         // the part type used for going around constraints
  part_t    start_trans; // for rubber treads, transition from arc to tangent
  part_t    end_trans;   // for rubber treads, transition from tangent to arc
} band_attrib_t;

extern band_attrib_t band_types[];
extern int n_band_types;
extern part_t band_constraints[];
extern int n_band_constraints;

void list_band_types(void);
void band_ini(       void);
int isbandtype(      char *type);
int isbandconstraint(char *type);

int
synth_band(
  char                *type,
  int                  n_constraints,
  LSL_band_constraint *constraints,
  int                  color,
  FILE                *output,
  int                  ghost,
  char                *group);

#ifdef _cplusplus
};
#endif

#endif
