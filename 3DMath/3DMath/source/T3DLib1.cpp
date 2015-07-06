//
//  T3DLib1.cpp
//  3DMath
//
//  Created by TTc on 15/3/7.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#include "T3DLib1.h"

#include <math.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdarg.h>
#include <ctype.h>

//GLOBALS
FILE    *fp_error  = NULL;      //general error file
char    error_filename[80];     //error file name









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


/**
 *  this function computers the distance from 0,0 to x,y with 3.5% error
 *  first computer the absoulute value of x,y
 *  @param x point x
 *  @param y point y
 *
 *  @return distance_2d
 */
int
Fast_Distance_2D(int x, int y) {
    x = abs(x);
    y = abs(y);
    int min = MIN(x,y);
    return (x+y-(min>>1)-(min>>2)+ (min>>4));
}

/**
 *  this function computers the distance from the origin to x,y,z
 *
 *  @param fx point x
 *  @param fy point y
 *  @param fz point z
 *
 *  @return distance_3D
 */
float
Fast_Distance_3D(float fx, float fy, float fz) {
    int temp; // used for swaping
    int x,y,z; //used for algorithm
    
    //make sure values are all postitive
    x = fabs(fx) * 1024;
    y = fabs(fy) * 1024;
    z = fabs(fz) * 1024;
    //sort values
    if (y < x) SWAP(x,y,temp);
    if (z < y) SWAP(y,z,temp);
    if (y < x) SWAP(x,y,temp);
    
    int dist = (z + 11 *(y >> 5) + (x >> 2));
    return (float)(dist >> 10);
}


int
Write_Error(char *string,...) {
    char buffer[1024];//working buffer
    
    va_list arglist;
    //make sure both the error file and string are valid
    if(!string || !fp_error)
        return (0);
    
    // print out the string using the variable number of arguments on stack
    va_start(arglist,string);
    vsprintf(buffer,string,arglist);
    va_end(arglist);
    fprintf(fp_error,buffer);
    
    // flush buffer incase the system bails
    fflush(fp_error);
    
    // return success
    return(1);

}









