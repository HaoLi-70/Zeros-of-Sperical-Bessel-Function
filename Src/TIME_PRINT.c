
#include "TIME_PRINT.h"

/**********************************************************************************************/

extern int Time_Print(void){
    
    /******************************************************************************************
     Purpose:
     Compute and print the running time.
     Record of revisions:
     30 Nov. 2019, Hao Li
     ******************************************************************************************/
    
    static int conts = 0;
    clock_t time_tmp;
    static clock_t time_begin = 0;
    
    if (conts == 0) {
        time_begin = clock();
        fprintf(stderr, "\n Time Calculating Is Initialized \n");
        conts++;
        
    }else{
        time_tmp = clock();
        long hours = (time_tmp-time_begin)/CLOCKS_PER_SEC/3600;
        
        long minutes = ((time_tmp-time_begin)/CLOCKS_PER_SEC-hours*3600)/60;
        
        double seconds = (time_tmp-time_begin)*1.0 \
            /CLOCKS_PER_SEC-hours*3600-minutes*60;
        
        fprintf(stderr, "\n Time Print Point %d \n",conts);
        
        if (hours > 0) {
            fprintf(stderr," Running time= %lu h %lu min %.2lf sec\n", \
                    hours,minutes,seconds);
        }else if(minutes > 0){
            fprintf(stderr," Running time= %lu min %.2lf sec\n",minutes,seconds);
        }else{
            fprintf(stderr," Running time= %.2lf sec\n",seconds);
        }
        conts++;
    }
    return conts-1;
}

/**********************************************************************************************/

