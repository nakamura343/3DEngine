//
//  T3DLib.c
//  3DMathLib
//
//  Created by TTc on 15/7/2.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#include "T3DLib4.h"
#include "T3DLib1.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>




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
    //提取角度的整数部分和小数部分 ,以便进行插值计算
    int  theta_int = (int)theta;
    int  theta_frac = theta - theta_int;
    
    return (sin_look[theta_int] +
            theta_frac * (sin_look[theta_int + 1] - sin_look[theta_int]));
}


float
Fast_Cos(float theta) {
    theta = fmodf(theta, 360);
    
    if(theta < 0)    theta += 360.0;
    
    int theta_int = (int)theta;
    float theta_frac = theta - theta_int;
    
    return (cos_look[theta_int] +
            theta_frac * (cos_look[theta_int +1] - cos_look[theta_int]));
}
/******************************************************************************/


/******************************************************************************/


//2D极坐标  x = r* sin(theta)
/*
 将一个用r和 theta表示的2D极坐标 点转换为一个(x,y)点,并将它存储在一个POINT2D结构中
 */
void
POLAR2D_To_POINT2D(POLAR2D_PTR polar, POINT2D_PTR rect) {
    rect->x = polar->r * cosf(polar->theta);
    rect->y = polar->r * sinf(polar->theta);
}
/*
 将一个2D极坐标转换为一个x,y坐标
 */
void
POLAR2D_To_RectXY(POLAR2D_PTR polar, float *x, float *y) {
    *x = polar->r * cosf(polar->theta);
    *y = polar->r * sinf(polar->theta);
}
/*
  将一个笛卡尔坐标系中的2D点(直角坐标系)转换为一个2D极坐标格式
 */
void
POINT2D_To_POLAR2D(POINT2D_PTR rect,POLAR2D_PTR polar) {
    polar->r     = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    polar->theta = atanf(rect->y / rect->x);
}
/*
 将)一个2D极坐标格式  转换为 r 和theta
 */
void
POINT2D_To_PolarRTh(POINT2D_PTR rect,float *r , float *theta) {
    *r     = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    *theta = atanf(rect->y / rect->x);
}


//柱面坐标 函数
/*
 将一个3D柱面坐标点转换为一个3D 直角坐标系点 ;将圆柱转换为矩形
 */
void
CYLINDRICAL3D_To_POINT3D(CYLINDRICAL3D_PTR cyl,POINT3D_PTR rect) {
    rect->x = cyl->r * cosf(cyl->theta);
    rect->y = cyl->r * sinf(cyl->theta);
    rect->z = cyl->z;
}

/*
 将一个3D柱面坐标点转换为 x, y , z 坐标;将圆柱转换为矩形
 */
void
CYLINDRICAL3D_To_RectXYZ(CYLINDRICAL3D_PTR cyl,float *x,float *y,float *z) {
    *x = cyl->r * cosf(cyl->theta);
    *y = cyl->r * sinf(cyl->theta);
    *z = cyl->z;
}
/*
 将一个3D 直角坐标系点 转换为 一个3D柱面坐标点 ;将圆柱转换为矩形
 */
void
POINT3D_To_CYLINDRICAL3D(POINT3D_PTR rect,CYLINDRICAL3D_PTR cyl) {
    cyl->r      = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    cyl->theta  = atanf(rect->y / rect->x);
    cyl->z      = rect->z;
}
/*
 将一个3D 直角坐标系点 转换为 3D柱面坐标 r, theta, z ;将圆柱转换为矩形
 */
void
POINT3D_To_CylindricalRThZ(POINT3D_PTR rect,float *r ,float *theta, float *z) {
    *r      = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    *theta  = atanf(rect->y / rect->x);
    *z      = rect->z;
}


//球面坐标 函数
/*
     将一个3D球面坐标点转换为一个3D直角坐标点 ; 将球面转化为矩形
 */
void
SPHERICAL3D_To_POINT3D(SPHERICAL3D_PTR sph,POINT3D_PTR rect) {
    float r;
    r       = sph->p * sinf(sph->phi);
    rect->z = sph->p * cosf(sph->phi);
    //使用r简化计算x,y
    rect->x = r * cosf(sph->theta);
    rect->y = r * sinf(sph->theta);
}
/*
 将一个3D球面坐标点转换为 x, y, z坐标 ; 将球面转化为矩形
 */
void
SPHERICAL3D_To_RectXYZ(SPHERICAL3D_PTR sph, float *x,float *y,float *z) {
    float r;
    r       = sph->p * sinf(sph->phi);
    *z = sph->p * cosf(sph->phi);
    //使用r简化计算x,y
    *x = r * cosf(sph->theta);
    *y = r * sinf(sph->theta);
}
/*
 将一个3D直角坐标点 转换为 一个3D球面坐标点 ; 将球面转化为矩形
 */
void
POINT3D_To_SPHERICAL3D(POINT3D_PTR rect,SPHERICAL3D_PTR sph) {
    sph->p = sqrtf((rect->x * rect->x) + (rect->y * rect->y) + (rect->z * rect->z));
    sph->theta = atanf(rect->y/ rect->x);
    
    float r = sph->p * sinf(sph->phi);
    sph->phi = asinf(r/sph->p);
}
/*
 将一个3D直角坐标点  转换为 一个3D球面坐标 p ,theta 和 phi; 将球面转化为矩形
 */
void
POINT3D_To_SphericalPThPh(POINT3D_PTR rect,float *p ,float *theta, float *phi) {
   *p = sqrtf((rect->x * rect->x) + (rect->y * rect->y) + (rect->z * rect->z));
    *theta = atanf(rect->y/ rect->x);
    
    float r = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    *phi    = asinf(r / (*p));
}

/******************************************************************************/













/******************************************************************************/
/******************************************************************************/
/******************************************************************************/




















