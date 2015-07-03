//
//  T3DLib.c
//  3DMathLib
//
//  Created by TTc on 15/7/2.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#include "T3DLib.h"
#include <math.h>




/******************************************************************************/

extern float cos_look[361];
extern float sin_look[361]; //储存0-360度的值 [0...360]

/*
  这个函数使用查找表sin_look[];但能够通过插值计算负角度和小数角度的正弦;因此与查找表相比,精度
  更高,但是速度低些.
 
 */
float
Fast_Sin(float theta) {
    theta = fmodf(theta,360); //将角度转换为0-359的值
    if(theta < 0) theta += 360.0;
    
    int  theta_int = (int)theta;
    int  theta_frac = theta - theta_int;
    
    return (sin_look[theta_int] +
            theta_frac * (sin_look[theta_int + 1] - sin_look[theta_int]));
}
/******************************************************************************/
