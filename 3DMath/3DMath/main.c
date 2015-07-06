//
//  main.c
//  3DMathLib
//
//  Created by TTc on 15/7/2.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "T3DLib.h"

void
test_MATRIX3X3(){
    
    MATRIX3X3 m = {1,2,3,4,5,6,7,8,9};
    m.M[2][1] = 100.0;
    int i , j;
    for (i  = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            printf("m[%d][%d]===>%.0f \n",i,j,m.M[i][j]);
        }
    }
}

void
test_Fast_sin(){
    float sin60 = Fast_Sin(60.0);
    printf("sin60==%.2f",sin60);
}

int
main(int argc, const char * argv[]) {
    // insert code here...
    printf("Hello, World!\n");
    
    //rtest_MATRIX3X3();
//    test_Fast_sin();
    float sin60 = Fast_Sin(60.0);
    printf("sin60==%.2f",sin60);
    
    return 0;
}
