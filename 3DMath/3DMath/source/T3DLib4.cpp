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


//2D向量函数
/*
   将指定向量va和vb相加,然后将结果存储在vsum中
 */
void
VECTOR2D_ADD(VECTOR2D_PTR va, VECTOR2D_PTR vb, VECTOR2D_PTR vsum) {
    // this function adds va+vb and return it in vsum
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
}

/*
  向量相加函数的堆栈版本,将结果返回到堆栈中
 */
VECTOR2D
VECTOR2D_ADD(VECTOR2D_PTR va , VECTOR2D_PTR vb) {
    //this function adds va + vb and returns the result on the stack. 
    VECTOR2D vsum ;
    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    return (vsum);
}
/*
 this function subtracets va-vb and return it in vdiff
 */
void
VECTOR2D_Sub(VECTOR2D_PTR va, VECTOR2D_PTR vb, VECTOR2D_PTR vdiff) {
    vdiff->x = va->x - vb->x;
    vdiff->y = va->y - vb->y;
}

/*
 this function subtracets va-vb and returns the result on the stack
 */
VECTOR2D
VECTOR2D_Sub(VECTOR2D_PTR va, VECTOR2D_PTR vb) {
    VECTOR2D vdiff;
    
    vdiff.x = va->x - vb->x;
    vdiff.y = va->y - vb->y;
    return (vdiff);
}
/*
 this function scales a vector by the constant k,
 leaves the original unchanged ,and returns the result in vscaled
 */
void
VECTOR2D_Scale(float k,VECTOR2D_PTR va,VECTOR2D_PTR vscaled) {
    //multiply each component by scaling factor
    vscaled->x = k *va->x;
    vscaled->y = k *va->y;
}
/*
 this function scales a vector by the constant k,
 by scaling the original
 */
void
VECTOR2D_Scale(float k,VECTOR2D_PTR va) {
    //multiply each component by scaling factor
    va->x *= k;
    va->y *= k;
}
/* 点积 
 computes the dot product between va and vb
*/
float
VECTOR2D_Dot(VECTOR2D_PTR va,VECTOR2D_PTR vb) {
    return (va->x * vb->x) + (va->y * vb->y);
}


/**
 *  computers the magnitude of a vector ,slow
 *
 *  @param va a vector
 *
 *  @return length of vector
 */
float
VECTOR2D_Length(VECTOR2D_PTR va) {
    return (sqrtf(va->x * va->x + va->y * va->y));
}
/**
 *  computers the magnitude of a vector using an approximation,fast
 *
 *  @param va a vector
 *
 *  @return length of vector
 */
float
VECTOR2D_Length_Fast(VECTOR2D_PTR va) {
    return (float)Fast_Distance_2D(va->x, va->y);
}
/**
 *  normalizes the sent vector in place
 *
 *  @param va vector
 */
void
VECTOR2D_Normalize(VECTOR2D_PTR va) {
    //computer length
    float length = sqrtf(va->x * va->x + va->y * va->y);
    //test for zero length vector
    //if found return zero vector
    if(length < EPSILON_E5) return;
    
    float length_inv = 1/length;
    va->x = va->x * length_inv;
    va->y = va->y * length_inv;
}

/**
 *  normalizes the sent vector and returns thr reuslt in vn
 *
 *  @param va sent vector
 *  @param vn result
 */
void
VECTOR2D_Normalize(VECTOR2D_PTR va,VECTOR2D_PTR vn) {
    VECTOR2D_ZERO(vn);
    
    // compute length
    float length = (float)sqrtf(va->x*va->x + va->y*va->y );
    
    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;
    
    float length_inv = 1/length;
    
    // compute normalized version of vector
    vn->x = va->x*length_inv;
    vn->y = va->y*length_inv;
}
/**
 * this function creates a vector from two vectors (or points)  in 3D space
 */
void
VECTOR2D_Build(VECTOR2D_PTR init,VECTOR2D_PTR term,VECTOR2D_PTR result) {
    result->x = term->x - init->x;
    result->y = term->y - init->y;
}
/**
 *  this function returns the cosine of the  angle between two vectors.
 *  Note,we could computer the actual angle ; many many times, in further calcs we will want ultimately
 *  compute cos of the angle, so why not just leave it!
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR2D_CosTh(VECTOR2D_PTR va, VECTOR2D_PTR vb) {
    return (VECTOR2D_Dot(va, vb)/ (VECTOR2D_Length(va) * VECTOR2D_Length(vb)));
}
/**
 *  this function prints out a VECTOR2D
 *
 *  @param va   <#va description#>
 *  @param name <#name description#>
 */
void
VECTOR2D_Print(VECTOR2D_PTR va, char *name) {
    Write_Error("\n%s=[",name);
    for (int index=0; index<2; index++)
        Write_Error("%f, ",va->M[index]);
    Write_Error("]");
}


/******************************************************************************/
//3D向量函数

/**
 *  this function adds va+vb and return it in vsum
 *
 *  @param va   <#va description#>
 *  @param vb   <#vb description#>
 *  @param vsum <#vsum description#>
 */
void
VECTOR3D_ADD(VECTOR3D_PTR va, VECTOR3D_PTR vb, VECTOR3D_PTR vsum) {
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
    vsum->z = va->z + vb->z;
}
/**
 *  this function adds va+vb and return the result on the stack.
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
VECTOR3D
VECTOR3D_ADD(VECTOR3D_PTR va , VECTOR3D_PTR vb) {
    VECTOR3D vsum;
    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    vsum.z = va->z + vb->z;
    return vsum;
}
/**
 *  this function subtracts va-vb and return it in vdiff the stack
 *
 *  @param va    <#va description#>
 *  @param vb    <#vb description#>
 *  @param vdiff <#vdiff description#>
 */
void
VECTOR3D_Sub(VECTOR3D_PTR va, VECTOR3D_PTR vb, VECTOR3D_PTR vdiff) {
    vdiff->x = va->x - vb->x;
    vdiff->y = va->y - vb->y;
    vdiff->z = va->z - vb->z;
}

/**
 *   this function subtracts va-vb and returns the result on  the stack
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
VECTOR3D
VECTOR3D_Sub(VECTOR3D_PTR va, VECTOR3D_PTR vb) {
    VECTOR3D vdiff;
    vdiff.x = va->x - vb->x;
    vdiff.y = va->y - vb->y;
    vdiff.z = va->z - vb->z;
    return vdiff;
}
/**
 *  this function scales a vector by the constant k,and modifies the original
 *
 *  @param k  <#k description#>
 *  @param va <#va description#>
 */
void
VECTOR3D_Scale(float k,VECTOR3D_PTR va) {
    va->x*=k;
    va->y*=k;
    va->z*=k;
}
/**
 *  this function scales a vector by the constant k,leaves the original 
 *  unchanged, and returns the result ; in vres as well as on the stack;
 *  @param k       <#k description#>
 *  @param va      <#va description#>
 *  @param vscaled <#vscaled description#>
 */
void
VECTOR3D_Scale(float k,VECTOR3D_PTR va,VECTOR3D_PTR vscaled) {
     vscaled->x = k * va->x;
     vscaled->y = k * va->y;
     vscaled->z = k * va->z;
}
/**
 * 点积 computes the dot product between va and vb
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR3D_Dot(VECTOR3D_PTR va,VECTOR3D_PTR vb) {
    return(va->x * vb->x) + (va->y * vb->y) + (va->z * vb->z);
}
/**    叉积
 *  this function computes the cross product between va and vb,
 *  and returns the vector that is perpendicular to each in vn
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *  @param vn <#vn description#>
 */
void
VECTOR3D_Cross(VECTOR3D_PTR va,VECTOR3D_PTR vb,VECTOR3D_PTR vn) {
    vn->x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn->y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn->z =  ( (va->x * vb->y) - (va->y * vb->x) );
}
/**
 *  this function computes the cross product between va and vb
 *  and returns the vector that is perpendicular to each
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
VECTOR3D
VECTOR3D_Cross(VECTOR3D_PTR va,VECTOR3D_PTR vb) {
    VECTOR3D vn;
    vn.x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn.y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn.z =  ( (va->x * vb->y) - (va->y * vb->x) );
    return vn;
}
/**
 *  computes the magnitude of a vector, slow
 *
 *  @param va <#va description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR3D_Length(VECTOR3D_PTR va) {
    return (float)sqrtf(va->x*va->x + va->y*va->y + va->z*va->z);
}
/**
 *  computes the magnitude of a vector using an approximation ; very fast
 *
 *  @param va <#va description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR3D_Length_Fast(VECTOR3D_PTR va) {
    return Fast_Distance_3D(va->x, va->y, va->z) ;
}
/**
 *  normalizes the sent vector in placew
 *
 *  @param va <#va description#>
 */
void
VECTOR3D_Normalize(VECTOR3D_PTR va) {
    float length = sqrtf(va->x*va->x + va->y*va->y + va->z*va->z);
    
    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;
    
    float length_inv = 1/length;
    
    // compute normalized version of vector
    va->x*=length_inv;
    va->y*=length_inv;
    va->z*=length_inv;
}
/**
 *   normalizes the sent vector and returns the result in vn
 *
 *  @param va <#va description#>
 *  @param vn <#vn description#>
 */
void
VECTOR3D_Normalize(VECTOR3D_PTR va,VECTOR3D_PTR vn) {
    VECTOR3D_ZERO(vn);

    // compute length
    float length = VECTOR3D_Length(va);
    
    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;
    
    float length_inv = 1.0/length;
    
    // compute normalized version of vector
    vn->x = va->x*length_inv;
    vn->y = va->y*length_inv;
    vn->z = va->z*length_inv;
}
/**
 this function creates a vector from two vectors (or points) in 3D space
 */
void
VECTOR3D_Build(VECTOR3D_PTR init ,VECTOR3D_PTR term , VECTOR3D_PTR result) {
    result->x = term->x - init->x;
    result->y = term->y - init->y;
    result->z = term->z - init->z;
}
/**
 *  this function returns the cosine of the angle between,
 *  two vectors. Note, we could compute the actual angle,
 *  many many times, in further calcs we will want ultimately
 *  compute cos of the angle, so why not just leave it
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR3D_CosTh(VECTOR3D_PTR va, VECTOR3D_PTR vb) {
    return(VECTOR3D_Dot(va,vb)/(VECTOR3D_Length(va)*VECTOR3D_Length(vb)));
}
/**
 *  this function prints out a VECTOR3D
 *
 *  @param va   <#va description#>
 *  @param name <#name description#>
 */
void
VECTOR3D_Print(VECTOR3D_PTR va, char *name) {
    Write_Error("\n%s=[",name);
    for (int index=0; index<3; index++)
        Write_Error("%f, ",va->M[index]);
    Write_Error("]");
}


/******************************************************************************/
/******************************************************************************/




















