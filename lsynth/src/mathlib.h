#ifndef MATHLIB
#define MATHLIB

#define PRECISION float

void        vectorcp(PRECISION dst[3],PRECISION src[3]);
void      vectoradd3(PRECISION dst[3],PRECISION lft[3], PRECISION rht[3]);
void       vectoradd(PRECISION lft[3],PRECISION rht[3]);
void      vectorsub3(PRECISION dst[3],PRECISION lft[3], PRECISION rht[3]);
void       vectorsub(PRECISION lft[3],PRECISION rht[3]);
PRECISION  vectorlen(PRECISION vect[3]);
void      vectorrot3(PRECISION res[3],PRECISION src[3],PRECISION rot[3][3]);
void       vectorrot(PRECISION loc[3],PRECISION m[3][3]);

void        matrixcp(PRECISION dst[3][3],PRECISION src[3][3]);
void       matrixadd(PRECISION dst[3][3],PRECISION src[3][3]);
void      matrixadd3(PRECISION res[3][3],PRECISION lft[3][3],PRECISION rht[3][3]);
void     matrixmult3(PRECISION res[3][3],PRECISION lft[3][3],PRECISION rht[3][3]);
void      matrixmult(PRECISION res[3][3],PRECISION src[3][3]);
void       matrixinv(PRECISION inv[3][3],PRECISION src[3][3]);
void       matrixneg(PRECISION neg[3][3],PRECISION src[3][3]);
int         matrixeq(PRECISION dst[3][3],PRECISION src[3][3]);

#endif

