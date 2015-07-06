//
//  T3DLib1.cpp
//  3DMath
//
//  Created by TTc on 15/7/6.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#include "T3DLib1.h"

#include <math.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

float cos_look[361];
float sin_look[361];
/*
 * create sin/cos lookup table
 */
void
Build_Sin_Cos_Tables(void) {
    //generate the tables 0 - 360 inclusive
    for (int ang = 0; ang <= 360; ang++) {
        //convert ang to radians
        float theta = (float)ang * PI / (float)180;
        //insert next entry into table
        cos_look[ang] = cos(theta);
        sin_look[ang] = sin(theta);
    }
}

