#ifndef ALLOC_H
#define ALLOC_H

#define farmalloc(n) malloc((n))
#define farrealloc(i,j) realloc((i),(j))
#define farcoreleft() (0)

#endif /* ALLOC_H */

