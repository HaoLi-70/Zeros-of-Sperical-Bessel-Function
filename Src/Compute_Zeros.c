
#include <stdio.h>
#include "ALLOCATION.h"
#include "READ_INPUT.h"
#include "SPECIAL_FUNCTIONS.h"
#include "TIME_PRINT.h"


#define Path_Input "./input.dat"

int main(int argc, const char * argv[]) {
    
    Time_Print();
    
    STRUCT_INPUT *Input = (STRUCT_INPUT *)malloc(sizeof(STRUCT_INPUT));

    char *inputpath = Path_Input;

    RDINPUT(Path_Input, Input);

    int Zero_L=60, Zero_N=60;
  
    Zeros_Output(Input->Path_Output, Zero_L, Zero_N, Input->BC, 1);
    
    free(Input);

    Time_Print();

    return 0;
}
