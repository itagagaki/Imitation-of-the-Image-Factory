#define NULSYM		0
#define RESSYM		1
#define CSGSYM		2
#define ATTSYM		3

#define NOOP			(-1)
#define BOE			0
#define	OP_EQUAL		1	/* =            */
#define	OP_EQUNION		2	/* +=           */
#define	OP_EQINTER		3	/* &=           */
#define	OP_EQDIFF		4	/* -=           */
#define	OP_UNION		5	/* +            */
#define	OP_INTER		6	/* &            */
#define	OP_NOT			7	/* ~            */
#define	OP_DIFF			8	/* -            */
#define SYMLPAR			9	/* (            */
#define SYMRPAR			10	/* )            */
#define	SYMEXPTERM		11	/* ;            */
#define ADDATR			12
#define SYMMOVE			13	/* move         */
#define SYMSCALE		14	/* scale        */
#define SYMROTATE		15	/* rotate       */
#define	SYMCOMMA		16	/* ,            */
#define SYMLMPAR		17	/* {            */
#define SYMRMPAR		18	/* }            */
#define SYMLBRACKET		19	/* [            */
#define SYMRBRACKET		20	/* ]            */
#define	SYMAXIS			21	/* ->		*/
#define	SYMOAXIS		22	/* =>		*/
#define SYMCUBE			23	/* cube         */
#define SYMRECTANGLE		24	/* rectangle    */
#define SYMPLANE		25	/* plane        */
#define SYMSPHERE		26	/* sphere       */
#define SYMELLIPSOID		27	/* ellipsoid    */
#define SYMHYPERBOLOID		28	/* hyperboloid  */
#define SYMCONE			29	/* cone         */
#define SYMPARABOLOID		30	/* paraboloid   */
#define SYMCYLINDER		31	/* cylinder     */
#define SYMSUR2ND		32	/* sur2nd	*/
#define	SYMATTRIBUTE		33	/* attrib       */
#define SYMAMBIENCE		34	/* ambience     */
#define SYMLIGHT		35	/* light        */
#define SYMCAMERA		36	/* camera       */
#define SYMPUT			37	/* put          */
#define SYMDIFFUSIBILITY	38	/* diffusibility */
#define SYMSPECULAR		39	/* specular     */
#define SYMTRANSPARENCY		40	/* transparency	*/
#define SYMHIGHLIGHT		41	/* highlight	*/
#define SYMREFRACT		42	/* refract      */
#define SYMPARALLEL		43	/* parallel     */
#define SYMINCLUDE		44	/* include      */
#define SYMTEXTURE		45	/* texture      */
#define SYMAT			46	/* at           */
#define	SYMLUMINOSITY		47	/* luminosity	*/
#define	SYMPOSITION		48	/* position	*/
#define	SYMDIRECTION		49	/* direction	*/
#define	SYMTARGET		50	/* target	*/
#define	SYMZOOM			51	/* zoom		*/
#define	SYMPOINT		52	/* point	*/
#define	SYMFPOINT		53	/* fpoint	*/
#define	SYMSPOT			54	/* spot		*/
#define	SYMSKY			55	/* sky		*/
#define	SYMX			56	/* x		*/
#define	SYMY			57	/* y		*/
#define	SYMZ			58	/* z		*/

#define	SYMPLUS  OP_UNION
#define	SYMMINUS OP_DIFF

#define OPCD_TOP	1
#define OPCD_LAST	11
#define	AFFIN_TOP	13
#define	AFFIN_LAST	15
#define PRCD_TOP	23
#define PRCD_LAST	32

union _symval {
   char code;
   OBJECT *objptr;
   ATTR *atptr;
};

typedef struct symtag {
   char *symword;
   int symtype;
   union _symval symval;
   struct symtag *lower_sym;
   struct symtag *upper_sym;
} SYMTAB; 

extern void init_sym(); 
extern void search_sym();
extern SYMTAB *def_symbol();
