/*
 * This file describes the interface to the LDRAW synthesizable parts library.
 * Kevin Clague
 */
#ifndef TUBE_H
#define TUBE_H

#ifdef _cplusplus
extern "C" {
#endif

typedef struct {
  char      type[256];
  char      descr[256];
  int       fill;         // FIXED or STRETCH
  PRECISION diameter;     // of cross section
  PRECISION stiffness;    // stiffness
  PRECISION twist;        // in degrees
  part_t    start;        // LDraw part for start of hose
  part_t    mid;          // LDraw part for mid section
  part_t    end;          // LDraw part for end of hose
  part_t    alt;          // LDraw part alternate for mid of hose
} hose_attrib_t;

extern hose_attrib_t hose_types[];
extern int n_hose_types;
extern part_t hose_constraints[];
extern int n_hose_constraints;

void list_hose_types(      void);
void list_hose_constraints(void);
void hose_ini(             void);
int ishosetype(      char *type);
int ishoseconstraint(char *type);

int
synth_hose(
  char   *type,
  int     n_constraints,
  part_t *constraints,
  int     ghost,
  char   *group,
  int     group_size,
  int     color,
  FILE   *output);
#ifdef _cplusplus
};
#endif

#endif
