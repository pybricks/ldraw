
#include "math.h"
#include "mathlib.h"

void
vectorcp(
  PRECISION dst[3],
  PRECISION src[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    dst[i] = src[i];
  }
}

void
vectoradd3(
  PRECISION dst[3],
  PRECISION lft[3],
  PRECISION rht[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    dst[i] = lft[i] + rht[i];
  }
}

void
vectoradd(
  PRECISION dst[3],
  PRECISION src[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    dst[i] += src[i];
  }
}

void
vectorsub3(
  PRECISION dst[3],
  PRECISION lft[3],
  PRECISION rht[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    dst[i] = lft[i] - rht[i];
  }
}

void
vectorsub(
  PRECISION dst[3],
  PRECISION src[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    dst[i] -= src[i];
  }
}

PRECISION
vectorlen(
  PRECISION vect[3])
{
  PRECISION len = 0;
  int i;

  for (i = 0; i < 3; i++) {
    len += vect[i]*vect[i];
  }

  return sqrt(len);
}

void
vectorrot3(
  PRECISION t[3],
  PRECISION r[3],
  PRECISION m[3][3])
{
#if 0
  t[0] = r[0]*m[0][0] + r[1]*m[1][0] + r[2]*m[2][0];
  t[1] = r[0]*m[0][1] + r[1]*m[1][1] + r[2]*m[2][1];
  t[2] = r[0]*m[0][2] + r[1]*m[1][2] + r[2]*m[2][2];
#else
  t[0] = r[0]*m[0][0] + r[1]*m[0][1] + r[2]*m[0][2];
  t[1] = r[0]*m[1][0] + r[1]*m[1][1] + r[2]*m[1][2];
  t[2] = r[0]*m[2][0] + r[1]*m[2][1] + r[2]*m[2][2];
#endif
}

void
vectorrot(
  PRECISION loc[3],
  PRECISION m[3][3])
{
  PRECISION t[3];

  vectorrot3(t,loc,m);

  loc[0] = t[0];
  loc[1] = t[1];
  loc[2] = t[2];
}

void
matrixcp(
  PRECISION dst[3][3],
  PRECISION src[3][3])
{
  int i,j;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dst[i][j] = src[i][j];
    }
  }
}

int
matrixeq(
  PRECISION dst[3][3],
  PRECISION src[3][3])
{
  int i,j,rc = 1;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      rc &= fabs(dst[i][j] - src[i][j]) < 0.0001;
    }
  }
  return rc;
}

void
matrixneg(
  PRECISION dst[3][3],
  PRECISION src[3][3])
{
  int i,j;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dst[i][j] = -src[i][j];
    }
  }
}

void
matrixadd3(
  PRECISION dst[3][3],
  PRECISION lft[3][3],
  PRECISION rht[3][3])
{
  int i,j;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dst[i][j] = lft[i][j] + rht[i][j];
    }
  }
}

void
matrixadd(
  PRECISION dst[3][3],
  PRECISION src[3][3])
{
  int i,j;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dst[i][j] += src[i][j];
    }
  }
}

void
matrixmult3(
  PRECISION res[3][3],
  PRECISION lft[3][3],
  PRECISION rht[3][3])
{
  int i,j,k;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      res[i][j] = 0.0;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        res[i][j] += lft[i][k] * rht[k][j];
      }
    }
  }
}

void
matrixmult(
  PRECISION res[3][3],
  PRECISION src[3][3])
{
  PRECISION t[3][3];

  matrixcp(t,res);
  matrixmult3(res,t,src);
}

void
matrixinv2(
  PRECISION inv[3][3],
  PRECISION src[3][3])
{
  int i;
  PRECISION tmp[3][3];

  for (i = 0; i < 3; i++) {
    int j;
    for (j = 0; j < 3; j++) {
      inv[i][j] = 0;
    }
    inv[i][i] = 1;
  }

  matrixcp(tmp,src);

  for (i = 0; i < 3; i++) {
    int       j;
    PRECISION l;

    l = tmp[i][i];

    if (l == 0) {
      break;
    } else {
      int a;
      for (j = 0; j <= 3-1; j++) {
        tmp[i][j]   /= l;
        inv[i][j] /= l;
      }
      for (a = 0; a < 3; a++) {
        if ((a-i) != 0) {
          PRECISION b;
          b = tmp[a][i];
          for (j = 0; j < 3; j++) {
            tmp[a][j]   -= b*tmp[i][j];
            inv[a][j] -= b*inv[i][j];
          }
        }
      }
    }
  }
}

void
matrixinv(
  PRECISION a[3][3],
  PRECISION src[3][3])
{
  int i,j,k;
  int p[3];
  PRECISION h,q,s,sup,pivot;

  matrixcp(a,src);

  for (k = 0; k < 3; k++) {
    sup = 0.0;
    p[k] = 0;
    for (i = k; i < 3; i++) {
      s = 0.0;
      for (j = k; j < 3; j++) {
        s += fabs(a[i][j]);
      }
      q = fabs(a[i][k])/s;
      if (sup < q) {
        sup = q;
        p[k] = i;
      }
    }
    if (sup == 0.0) {
      return;
    }
    if (p[k] != k) {
      for (j = 0; j < 3; j++) {
        h = a[k][j];
        a[k][j] = a[p[k]][j];
        a[p[k]][j] = h;
      }
    }
    pivot = a[k][k];
    for (j = 0; j < 3; j++) {
      if (j != k) {
        a[k][j] = -a[k][j] / pivot;
        for (i = 0; i < 3; i++) {
          if (i != k) {
            a[i][j] += a[i][k] * a[k][j];
          }
        }
      }
    }
    for (i = 0; i < 3; i++) {
      a[i][k] /= pivot;
    }
    a[k][k] = 1/pivot;
  }
  for (k = 3 - 1; k >= 0; k--) {
    if (p[k] != k) {
      for (i = 0; i < 3; i++) {
        h = a[i][k];
        a[i][k] = a[i][p[k]];
        a[i][p[k]] = h;
      }
    }
  }
}
