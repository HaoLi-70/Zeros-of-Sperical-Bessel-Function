
#ifndef RD_INPUT_h
#define RD_INPUT_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "STR.h"
#include "ALLOCATION.h"

enum keywordtype {KEYWORD_REQUIRED, KEYWORD_DEFAULT, KEYWORD_OPTIONAL};

typedef struct Struct_Keywords{
    char keyword[Key_Length];
    char line[Max_Line_Length];
    bool Set, Required;
}STRUCT_KEYS;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Input{

    //N and L
    int Zero_L, Zero_N;

    //0 zeros, 1 zero derivatives
    int BC;

    //Path to the output
    char Path_Output[Max_Line_Length];
    
}STRUCT_INPUT;

/*--------------------------------------------------------------------------------*/

extern void Keywords_Conversion(STRUCT_KEYS Keywords[], STRUCT_INPUT *Input);

extern void RDINPUT(char Filename[], STRUCT_INPUT *Input);

/*--------------------------------------------------------------------------------*/

#endif /* RD_INPUT_h */
