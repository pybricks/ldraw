#include <ctype.h>

int strncasecmp(
  const char *s1,
  const char *s2,
  int   n)
{
  unsigned c1, c2;

  if (n == 0) {
    return 0;
  }

  do {
    c1 = tolower(*s1++);
    c2 = tolower(*s2++);

    if (c1 != c2)
      return c1 - c2;
    if (c1 == 0)
      break;
  } while (--n != 0);
  return 0;
}

int strcasecmp(
  const char *s1,
  const char *s2)
{
  unsigned c1, c2;

  do {
    c1 = tolower(*s1++);
    c2 = tolower(*s2++);

    if (c1 != c2)
      return c1 - c2;
    if (c1 == 0)
      break;
  } while(1);

  return 0;
}
