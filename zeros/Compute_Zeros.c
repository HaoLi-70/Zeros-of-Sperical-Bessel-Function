
#include <stdio.h>
#include "SPECIAL_FUNCTIONS.h"
#include "TIME_PRINT.h"
#include "ALLOCATION.h"

#define Path_Zeros "/Users/Hao/Desktop/Zeros.dat"

#define Path_Derivative_Zeros "/Users/Hao/Desktop/Derivative_Zeros.dat"

int main(int argc, const char * argv[]) {
    
    Time_Print();
    
    char *filename;
    
    int Zero_L=60, Zero_N=60;
  
    filename=Path_Zeros;
    Zeros_Output(filename, Zero_L, Zero_N, 0, 1);
    filename=Path_Derivative_Zeros;
    Zeros_Output(filename, Zero_L, Zero_N, 1, 1);
    
    Time_Print();
    return 0;
}
