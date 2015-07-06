//
//  test_mathlib.cpp
//  3DMath
//
//  Created by TTc on 15/7/6.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#include "test_mathlib.h"

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
    Build_Sin_Cos_Tables();
    float sin60 = Fast_Sin(60);
    printf("sin60==%.2f",sin60);
}
void
test_Fast_cos(){
    Build_Sin_Cos_Tables();
    float cos60 = Fast_Cos(60);
    printf("cos60==%.2f",cos60);
}
