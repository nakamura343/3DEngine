//
//  T3DLib1.h
//  3DMath
//
//  Created by TTc on 15/3/7.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#ifndef ___DMath__T3DLib1__
#define ___DMath__T3DLib1__

#include <stdio.h>
//与Pi相关的常量
#define     PI          ((float)3.141592654f)
#define     PI2         ((float)6.283185307f)
#define     PI_DIV_2    ((float)1.570796327f)
#define     PI_DIV_4    ((float)0.785398163f)
#define     PI_INV      ((float)0.318309886f)


// used to compute the min and max of two expresions
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

//used for swapping algorithm
#define SWAP(a,b,t) {t=a; a=b; b=t;}

//some math macros
#define DEG_TO_RAD(ang) ((ang)*PI/180.0)
#define RAD_TO_DEG(rads) ((rads)*180.0/PI)

#define RAND_RANGE(x,y) ( (x) + (rand()%((y)-(x)+1)))

void Build_Sin_Cos_Tables(void);
// math functions
int Fast_Distance_2D(int x, int y);

float Fast_Distance_3D(float fx, float fy, float fz);

//error functions
int Write_Error(char *string,...);


#endif /* defined(___DMath__T3DLib1__) */
