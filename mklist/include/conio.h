#ifndef CONIO_H
#define CONIO_H

#define getch() getc(stdin)
  
int stricmp(const char *str1, const char *str2) 
{
  const unsigned char * ptr1 = (unsigned char *) str1;
  const unsigned char * ptr2 = (unsigned char *) str2;
  unsigned char c1, c2;
	
  while ((c1 = toupper(*ptr1++)) == (c2 = toupper(*ptr2++)))
  {
    if (!c1)
    {
      return(0); // end of both strings reached, so they match
    }
  }
  // first non-matching char was reached, including possibly 0 on one or the other string
  return(c1 - c2);
}

#endif /* CONIO_H */
