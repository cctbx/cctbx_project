/*    suitenscrt.h  see PKINSCRT.c    */

/*suitescrt.c defines SCRATCH as nothing ("") so really does the declarations*/
/* if SCRATCH not defined, then it is defined as "extern"  so only referenced*/

#ifdef  SCRATCH
#undef  SCRATCH
#define SCRATCH
#else
#define SCRATCH extern
#endif

SCRATCH  char alertstr[256];

typedef struct{
    char* begin;
    char* end;
    char* next;
    char* cursor;
} textblock;

SCRATCH textblock mainscratch;

/*prototypes*/
void   inittextblock(textblock*);
void   rewindtextblock(textblock*);
void   buildtextblock(textblock*, size_t);
void   disposetextblock(textblock*);
void   putonetextblockline(textblock*, char*);
void   getonetextblockline(textblock*, char*);


