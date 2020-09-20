#ifndef casecmp_h
#define casecmp_h

#ifdef WIN32
#ifdef __MINGW_H
  // These are in string.h in mingw version of gcc. 
  // Others may need to include <strings.h> but it usually conflicts with <string.h>
  // Maybe we should move the #include <string.h> from lsynthcp.h to here.
#else
int strncasecmp(
  const char *s1,
  const char *s2,
  int   n);

int strcasecmp(
  const char *s1,
  const char *s2);
#endif
#endif

#endif
