#include <stdio.h>

int putd(value, fp)
double value;
FILE *fp;
{
  register char *ptr, i;

  ptr = (char *)&value;
  for (i = 0; i < sizeof(double); i++)
    if (putc(*ptr++, fp) == EOF) return EOF;

  return 0;
}

