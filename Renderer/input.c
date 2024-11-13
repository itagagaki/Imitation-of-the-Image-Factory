#include <stdio.h>

double getd(stream)
FILE *stream;
{
  double value;
  register char *ptr, i;

  ptr = (char *)&value;
  for (i = 0; i < sizeof(double); i++) *ptr++ = getc(stream);

  return value;
}
