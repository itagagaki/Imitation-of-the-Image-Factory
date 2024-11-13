#define MAXSYMLEN 31

extern int lexi();
extern int fget_token();
extern void check_token();
extern bool finput();

extern char string[];
extern char word[];
extern int toktype;
extern int opcode;
extern SYMTAB *sym;
extern int ch;
extern FILE *rfp;
extern unsigned long *inpline;
