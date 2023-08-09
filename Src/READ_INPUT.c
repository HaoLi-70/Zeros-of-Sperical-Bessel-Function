
#include "READ_INPUT.h"

/*--------------------------------------------------------------------------------*/

  /*######################################################################
   
    revision log:
        9 Aug. 2023.
   
   ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern void Keywords_Conversion(STRUCT_KEYS Keywords[], STRUCT_INPUT *Input){
  
    /*######################################################################
      Purpose:
        convert the keywords to configuration.
      Record of revisions:
        9 Aug. 2023.
      Input parameters:
        Keywords, a structure with the keywords.
      Output parameters:
        Input, a structure with the input information.
     ######################################################################*/
        
    int indx;

    indx = 0;
    Input->BC = atoi(Keywords[indx].line);
    if(Input->BC==0){
      fprintf(stderr, "\n compute location of zero values \n");
    }else{
      fprintf(stderr, "\n compute location of zero derivatives \n");
    }

    indx = 1;
    Input->Zero_L = atoi(Keywords[indx].line);
    fprintf(stderr, "\n L = %d \n", Input->Zero_L);

    indx = 2;
    Input->Zero_N = atoi(Keywords[indx].line);
    fprintf(stderr, "\n N = %d \n", Input->Zero_N);

    indx = 3;
    String_Copy(Input->Path_Output, Keywords[indx].line, \
        strlen(Keywords[indx].line), true);
    fprintf(stderr, "\n output path : %s \n", Input->Path_Output);
    
    return;
}

/*--------------------------------------------------------------------------------*/

extern void RDINPUT(char Filename[], STRUCT_INPUT *Input){
  
    /*######################################################################
      Purpose:
        Read the input file for zero position calculation.
      Record of revisions:
        9 Aug. 2023.
      Input parameters:
        Filename[], the input file.
      Output parameters:
        Input, a structure saved the input information.
     ######################################################################*/
    
    const char *routine_name = "RDINPUT";
    
    FILE *fa;
    fa = fopen(Filename, "r");
        
    int Num_Keywords = 3;

    STRUCT_KEYS Keywords[] ={
      {"boundary_condition", "", false, true}, //0
      {"Zero_L", "60", false, false}, //1
      {"Zero_N", "60", false, false},  //2
      {"Path_Output", "./zeros", false, false},  //3
    };
    Num_Keywords = sizeof(Keywords)/sizeof(STRUCT_KEYS);
    
    char lines[Max_Line_Length], key[Key_Length], *p;
    int len, i;
    long len_tot;
    bool neglect = true;
    
    while (Read_line(lines, fa) > 0 ) {
      len_tot = strlen(lines);
      len = Indx_Char(lines, '=', 1);
      if (len > 1){
        len_tot -= len;
        p = lines+len;
        String_Copy(key, lines, len-1, true);
        Trim(key, 3);
        
        for (i=0; i<Num_Keywords; i++) {
          if (strcmp(key,Keywords[i].keyword)==0){
            if (Keywords[i].Set == true){
              fprintf(stderr, "warning: read keyword %s \
                  more than once \n", key);
            }
            String_Copy(Keywords[i].line, p, len_tot, true);
            Keywords[i].Set = true;

            fprintf(stderr," %d %s %s \n",i, Keywords[i].keyword, \
                Keywords[i].line);
            
            neglect = false;
            break;
          }
        }
        
        if (neglect) {
          fprintf(stderr, "Warning : Neglect key words %s \n",key);
          fprintf(stderr, "\n The neglected line is %s \n",lines);
        }
      }
    }
    
    fclose(fa);

    Keywords_Conversion(Keywords, Input);

    for (i=0; i<Num_Keywords; i++) {
      if (Keywords[i].Required == true && Keywords[i].Set == false){
        fprintf(stderr, "Error : do not find keyword %s \n", \
            Keywords[i].keyword);
        exit(1);
      }
    }

  return;
}

/*--------------------------------------------------------------------------------*/


