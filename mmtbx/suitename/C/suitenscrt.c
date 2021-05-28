/*    suitenscrt.c  cribbed from PKINSCRT.c  070210  */

#include "suitename.h"
#define SCRATCH
    /*that is, SCRATCH is blank, so suitenscrt.h will do initial declarations*/
#include "suitenscrt.h"
#undef  SCRATCH

/****inittextblock()**********************************************************/
void inittextblock(textblock* tb)
{
    tb->begin = tb->end = tb->next = tb->cursor = NULL;
}
/*___inittextblock()_________________________________________________________*/

/****rewindtextblock()********************************************************/
void rewindtextblock(textblock* tb)
{
    tb->cursor = tb->begin; /*just move cursor back to beginning*/
}
/*___rewindtextblock()_______________________________________________________*/

/*3456789_123456789_123456789_123456789_123456789_123456789_123456789_12345678*/
/****disposetextblock()*******************************************************/
void disposetextblock(textblock* tb)
{
    if (tb->begin) 
    {
        free(tb->begin);   /*actually releases any allocated space*/
        inittextblock(tb); /*just sets member pointers==0*/
    }
}
/*___disposetextblock()______________________________________________________*/
/*3456789_123456789_123456789_123456789_123456789_123456789_123456789_12345678*/
/****putonetextblockline()****************************************************/
void putonetextblockline(textblock* tb, char* thestring)
{
    static char stringtemp[256];

    /*store directly in giant character array, textblock */
    if( (unsigned)(tb->end - tb->next) < strlen(thestring)+8 ) /*give some slack*/
    {/*build more space*/
        strcpy(stringtemp,thestring);/*leaving this routine forgets contents*/
        buildtextblock(tb, 10000); /*PKINSCRT.c*/ /*try for (n) more */
        strcpy(thestring,stringtemp);
    }
    if( (unsigned)(tb->end - tb->next) > strlen(thestring) ) /*insurance */
    {
        strcpy(tb->next,thestring);
        tb->next += strlen(tb->next) + 1;
    }
}
/*___putonetextblockline()___________________________________________________*/

/****getonetextblockline()****************************************************/
void getonetextblockline(textblock* tb, char* thestring)
{
    /*get one character string from textblock, giant character array*/
    if(tb->cursor < tb->next)  /*next available open space*/
    {
        strcpy(thestring,tb->cursor);
        tb->cursor += strlen(tb->cursor) + 1;
    }
    else thestring[0] = '\0'; /*NULL string signel for end of data*/
    return;
}
/*___getonetextblockline()___________________________________________________*/

/****buildtextblock()*********************************************************/
void buildtextblock(textblock* tb, size_t more)
{
    textblock tbtemp;
        /*members: char *begin, *end, *next, *cursor */    
    size_t length,nextoffset;
    int  ifail, recycle, moretry;
    
    tbtemp.begin = NULL;
    tbtemp.end = NULL;
    tbtemp.next = NULL;
    tbtemp.cursor = NULL;

    moretry = more;
    recycle = 1;
    while(recycle==1)
    {/*recycle allocation trials*/
        ifail = 0;
        length = tb->end  - tb->begin;
        nextoffset = tb->next - tb->begin;
        if(tb->begin)
        {/*exists: realloc more space*/
           length = tb->end - tb->begin;
           tbtemp.begin=(char *)realloc(tb->begin,sizeof(char)*(length+more));
           if(tbtemp.begin==NULL) ifail = 1;
        }
        else
        {/*not previously built, malloc initial "more" space*/
            length = 0;
            nextoffset = 0;
            tbtemp.begin = (char *)malloc(sizeof(char)*more);
            if(tbtemp.begin==NULL) ifail = 1;
        }
        if(ifail)
        {/*failed to do a reallocation, reduce request and try again*/
            more = more/2;
            if(more > 3) recycle = 1;
            else         recycle = 0;
        }
        else recycle = 0;
    }/*recycle allocation trials*/

    if(ifail==0)  
    {
        tb->begin = tbtemp.begin; /*reassign pointer*/
        tb->next = tb->begin + nextoffset; /*reset*/
        tb->cursor = tb->next; /*current position*/
        tb->end = tb->begin + sizeof(char)*(length+more);
    }
    if(Ltest||ifail)
    {
        if(ifail)
            sprintf(alertstr,CRLF"+%d text block reallocation failed"CRLF
                             ", remains: %ld"
                    ,moretry,(tb->end - tb->begin));
        else
            sprintf(alertstr,CRLF"text block allocation now == %ld"CRLF
                   ,(tb->end - tb->begin));
        fprintf(stderr,"%s",alertstr);
    }
}
/*___buildtextblock()________________________________________________________*/

