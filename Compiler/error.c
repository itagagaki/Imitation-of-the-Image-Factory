#include <stdio.h>
#include "config.h"
#include "error.h"

extern char *inpfile, *outfile;
extern unsigned long *inpline;
extern FILE *wfp;

static char *errmsg[] = {
  "",
  "bad syntax",				/*  1: BAD_SYNTAX	*/
  "out of primitive table",		/*  2: TOO_SHAPE	*/
  "invalid primitive name",		/*  3: IN_PRIMNAME	*/
  "too many attributes",		/*  4: TOO_ATTR		*/
  "too many objects",			/*  5: TOO_OBJECT	*/
  "too many lights",			/*  6: TOO_LIGHT	*/
  "multiple camera descriptions",	/*  7: MULT_CAM		*/
  "not constant",			/*  8: NOT_CONST	*/
  "bad constant",			/*  9: BAD_CONST	*/
  "unexpected EOF in comment",		/* 10: EOF_IN_COMMENT	*/
  "bad character,",			/* 11: BAD_CHR		*/
  "invalid operator",			/* 12: INVALOP		*/
  "multiple camera position",		/* 13: MULT_CPOS	*/
  "multiple camera direction",		/* 14: MULT_CDIR	*/
  "multiple camera rotate",		/* 15: MULT_CROT	*/
  "multiple camera zoom",		/* 16: MULT_CZOOM	*/
  "multiple camera iris",		/* 17: MULT_CIRIS	*/
  "missing camera position",		/* 18: NO_CAMPOS	*/
  "multiple sky",			/* 19: MULT_SKY		*/
  "bad axis of rotation",		/* 20: BAD_AXIS		*/
  "invalid omossion",			/* 21: INV_OMIT		*/
  "can't open input",			/* 22: CANT_INPUT	*/
  "can't create",			/* 23: CANT_CREATE	*/
  "string too long",			/* 24: STR_TOOLONG	*/
  "file name too long",			/* 25: STR_TOOLONG	*/
  "string needs terminating \"",	/* 26: UNTERM_STRING	*/
};

static void cancel()
{
  if (wfp) fclose(wfp);
  unlink(outfile);
}

void int_error(msg)
char *msg;
{
  fprintf(stderr, "internal error - %s\n", msg);
  cancel();
  exit(1);
}

void error(level, code, arg)
int level;
int code;
char *arg;
{
  fprintf(stderr, "%s", inpfile);
  if (inpline) fprintf(stderr, " %u", *inpline);
  fprintf(stderr, ": %s: %s", level ? "ERROR" : "warning", errmsg[code]);
  if (arg) fprintf(stderr, " %s", arg);
  fputc('\n', stderr);
  if (level) {
    cancel();
    exit(1);
  }
}
