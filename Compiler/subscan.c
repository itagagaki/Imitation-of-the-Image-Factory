#include <stdio.h>
#include <ctype.h>
#include "config.h"
#include "vector.h"
#include "object.h"
#include "symbol.h"
#include "token.h"
#include "subscan.h"
#include "error.h"

#define MAXSTRLEN 255

FILE *rfp;
unsigned long *inpline;
int ch;
int toktype;
int opcode;
SYMTAB *sym;
char string[( MAXSTRLEN + 1 )];
char word[( MAXSYMLEN + 1 )];
double value;

int getch()
{
  ch = getc(rfp);

#if 0
  if (ch != EOF) putchar(ch);
#endif

  if (ch == '\n') ++(*inpline);
  return ch;
}

/*
<string> ::= "{<CHARACTER>}"
*/
static void read_string()
{
  int count;
  char *ptr;
  bool flag = TRUE;

  ptr = string;
  count = MAXSTRLEN;

  while (ch != '"') {
    if (ch == EOF || ch == '\n' || ch == '\r') error(1, UNTERM_STRING, NULL);

    if (count > 0) {
      *ptr++ = ch;
      --count;
    }
    else if (flag) {
      error(0, STR_TOOLONG, NULL);
      flag = FALSE;
    }

    getch();
  }

  getch();
  *ptr = '\0';
  toktype = STRTOK;
}

/*
<word>  ::= <alpha>{<alpha>|<num>}
<alpha> ::= A|B|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z
	   |a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r|s|t|u|v|w|x|y|z
	   |_
<num>   ::= 0|1|2|3|4|5|6|7|8|9
*/
static void read_word()
{
  int count;
  char *ptr;

  ptr = word;
  count = MAXSYMLEN;
  do {
    if (count > 0) {
      *ptr++ = ch;
      --count;
    }
    getch();
  } while (isalnum(ch) || ch == '_');

  *ptr = '\0';
  toktype = SYMTOK;
  search_sym(word);
}

/*
<number> ::= [<s>](<d>{<d>}[.{<d>}])|(.<d>{<d>})[e|E[<s>]{<d>}]
<s>      ::= +|-
<d>      ::= 0|1|2|3|4|5|6|7|8|9
*/
static void read_number(first_ch)
int first_ch;
{
  int sign = 1;
  double val = 0, power = 1;

  if (first_ch == '+' || first_ch == '-')
    sign = (first_ch == '+' ? 1 : -1);

  if (first_ch == '.')
    if (isdigit(ch)) {
      do {
	power *= 10;
	val = 10 * val + ch - '0';
	getch();
      } while (isdigit(ch));
    }
    else {
      toktype = BADTOK;
      return;
    }
  else {
    do {
      val = 10 * val + ch - '0';
      getch();
    } while (isdigit(ch));

    if (ch == '.') {
      while (isdigit(getch())) {
	power *= 10;
	val = 10 * val + ch - '0';
      }
    }
  }

  if (ch == 'e' || ch == 'E') {
    unsigned int i;
    int esign = 1;
    unsigned int eval = 0;
    double epower = 1.0;

    getch();
    if (ch == '+' || ch == '-') {
      esign = (ch == '+' ? 1 : -1);
      getch();
    }

    while (isdigit(ch)) eval = 10 * eval + ch - '0';

    for (i = 0; i < eval; i++) epower *= 10.0;

    if (esign == 1) value = sign * val / power * epower;
    else value = sign * val / (power * epower);
  }
  else value = sign * val / power;

  toktype = NUMTOK;
}

int lexi()
{
  int ch1;

re_lexi:
  while (isspace(ch)) getch();

  if (isalpha(ch) || ch == '_') read_word();
  else if (isdigit(ch)) read_number(NULL);
  else {
    ch1 = ch;
    getch();
    switch (ch1)
      {
      case '/':
	if (ch == '*') {
	  getch();
	  do {
	    if (ch == EOF) error(1, EOF_IN_COMMENT, NULL);
	    ch1 = ch;
	    getch();
	  } while (ch1 != '*' || ch != '/');

	  getch();
	  goto re_lexi;
	}
	else if (ch == '/') {
	  while (ch != EOF && ch != '\n' && ch != '\r') getch();
	  getch();
	  goto re_lexi;
	}
	else toktype = BADTOK;
	break;

      case '"':
	read_string();
	break;

      case '.':
	read_number(ch1);
	break;

      case '=':
	toktype = OP_TOK;
	if (ch == '>') {
	  opcode = SYMOAXIS;
	  getch();
	}
	else opcode = OP_EQUAL;
	break;

      case '+':
	if (isdigit(ch) || ch == '.') read_number(ch1);
	else {
	  toktype = OP_TOK;
	  if (ch == '=') {
	    opcode = OP_EQUNION;
	    getch();
	  }
	  else opcode = OP_UNION;
	}
	break;

      case '-':
	if (isdigit(ch) || ch == '.') read_number(ch1);
	else {
	  toktype = OP_TOK;
	  if (ch == '=') {
	    opcode = OP_EQDIFF;
	    getch();
	  }
	  else if (ch == '>') {
	    opcode = SYMAXIS;
	    getch();
	  }
	  else opcode = OP_DIFF;
	}
	break;

      case '&':
	toktype = OP_TOK;
	if (ch == '=') {
	  opcode = OP_EQINTER;
	  getch();
	}
	else opcode = OP_INTER;
	break;

      case '!':
	toktype = OP_TOK;
	opcode = OP_NOT;
	break;

      case '(':
	toktype = OP_TOK;
	opcode = SYMLPAR;
	break;

      case ')':
	toktype = OP_TOK;
	opcode = SYMRPAR;
	break;

      case '{':
	toktype = OP_TOK;
	opcode = SYMLMPAR;
	break;

      case '}':
	toktype = OP_TOK;
	opcode = SYMRMPAR;
	break;

      case '[':
	toktype = OP_TOK;
	opcode = SYMLBRACKET;
	break;

      case ']':
	toktype = OP_TOK;
	opcode = SYMRBRACKET;
	break;

      case ';':
	toktype = OP_TOK;
	opcode = SYMEXPTERM;
	break;

      case ',':
	toktype = OP_TOK;
	opcode = SYMCOMMA;
	break;

      case EOF:
	toktype = NULL;
	break;

      default:
	{
	  char temp[14];

	  sprintf(temp, " %c, hex 0X%02x", isprint(ch) ? ch : ' ', ch);
	  error(1, BAD_CHR, temp);
	}
      }
  }

  return toktype;
}

int fget_token()
{
  switch (lexi())
    {
    case SYMTOK:
      if (sym == NULL || sym->symtype != RESSYM) error(1, BAD_SYNTAX, NULL);
      return sym->symval.code;

    case OP_TOK:
      return opcode;
    }
}

void check_token(code)
int code;
{
  if (fget_token() != code) error(1, BAD_SYNTAX, NULL);
}

/* input numeric double precision value */
bool finput(val, token, defau_sw)
double *val;
int token;
bool defau_sw;
{
  switch (lexi())
    {
    case NUMTOK:
      *val = value;
      check_token(token);
      return FALSE;

    case OP_TOK:
      if (opcode != token) error(1, BAD_SYNTAX, NULL);
      if (defau_sw) return TRUE;
      else error(1, INV_OMIT, NULL);

    case SYMTOK:
      if (sym == NULL) error(1, NOT_CONST, word);
      if (sym->symtype != RESSYM) error(1, NOT_CONST, word);
      if (sym->symval.code != token) error(1, NOT_CONST, word);
      if (defau_sw) return TRUE;
      else error(1, INV_OMIT, NULL);

    default:
      error(1, BAD_SYNTAX, NULL);
    }
}

/* end */
