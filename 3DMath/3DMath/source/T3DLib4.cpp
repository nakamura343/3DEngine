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

/**
 *  this function adds va+vb and return it in vsum
 *
 *  @param va   <#va description#>
 *  @param vb   <#vb description#>
 *  @param vsum <#vsum description#>
 */
void
VECTOR4D_ADD(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum) {
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
    vsum->z = va->z + vb->z;
    vsum->w = 1;
}
/**
 *  this function adds va+vb and returns the result on  the stack
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
VECTOR4D
VECTOR4D_ADD(VECTOR4D_PTR va , VECTOR4D_PTR vb) {
    VECTOR4D vsum;
    
    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    vsum.z = va->z + vb->z;
    vsum.w = 1;
    return vsum;
}
/**
 *  this function subtracts va-vb and return it in vdiff
 *
 *  @param va    <#va description#>
 *  @param vb    <#vb description#>
 *  @param vdiff <#vdiff description#>
 */
void
VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vdiff) {
    vdiff->x = va->x - vb->x;
    vdiff->y = va->y - vb->y;
    vdiff->z = va->z - vb->z;
    vdiff->w = 1;
}
/**
 *  this function subtracts va-vb and returns the result on the stack
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
VECTOR4D
VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb) {
    VECTOR4D vdiff;
    
    vdiff.x = va->x - vb->x;
    vdiff.y = va->y - vb->y;
    vdiff.z = va->z - vb->z;
    vdiff.w = 1;
    return vdiff;
}
/**
 *  this function scales a vector by the constant k,in place , note w is left unchanged
 *
 *  @param k  <#k description#>
 *  @param va <#va description#>
 */
void
VECTOR4D_Scale(float k,VECTOR4D_PTR va) {
    // multiply each component by scaling factor
    va->x*=k;
    va->y*=k;
    va->z*=k;
    va->w = 1;
}
/**
 *  this function scales a vector by the constant k,
 *  leaves the original unchanged, and returns the result
 *  in vres as well as on the stack
 *  @param k       <#k description#>
 *  @param va      <#va description#>
 *  @param vscaled <#vscaled description#>
 */
void
VECTOR4D_Scale(float k,VECTOR4D_PTR va,VECTOR4D_PTR vscaled) {
    // multiply each component by scaling factor
    vscaled->x = k*va->x;
    vscaled->y = k*va->y;
    vscaled->z = k*va->z;
    vscaled->w = 1;
}
/**  点积
 *  computes the dot product between va and vb
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR4D_Dot(VECTOR4D_PTR va,VECTOR4D_PTR vb) {
    return (va->x * vb->x) + (va->y * vb->y) + (va->z * vb->z) ;
}
/**  叉积
 *  this function computes the cross product between va and vb
 *  and returns the vector that is perpendicular to each in vn
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *  @param vn <#vn description#>
 */
void
VECTOR4D_Cross(VECTOR4D_PTR va,VECTOR4D_PTR vb,VECTOR4D_PTR vn) {
    vn->x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn->y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn->z =  ( (va->x * vb->y) - (va->y * vb->x) );
    vn->w = 1;
}
/**
 *  this function computes the cross product between va and vb
 *  and returns the vector that is perpendicular to each
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
VECTOR4D
VECTOR4D_Cross(VECTOR4D_PTR va,VECTOR4D_PTR vb) {
    VECTOR4D vn;
    
    vn.x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn.y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn.z =  ( (va->x * vb->y) - (va->y * vb->x) );
    vn.w = 1;
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
VECTOR4D_Length(VECTOR4D_PTR va) {
    return sqrtf(va->x*va->x + va->y*va->y + va->z*va->z) ;
}
/**
 *  computes the magnitude of a vector using an approximation ;very fast
 *
 *  @param va <#va description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR4D_Length_Fast(VECTOR4D_PTR va) {
    return Fast_Distance_3D(va->x, va->y, va->z) ;
}
/**
 *  normalizes the sent vector and returns the result
 *
 *  @param va <#va description#>
 */
void
VECTOR4D_Normalize(VECTOR4D_PTR va) {
    // compute length
    float length = sqrtf(va->x*va->x + va->y*va->y + va->z*va->z);
    
    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;
    
    float length_inv = 1.0/length;
    
    // compute normalized version of vector
    va->x*=length_inv;
    va->y*=length_inv;
    va->z*=length_inv;
    va->w = 1;
}
/**
 *  normalizes the sent vector and returns the result in vn
 *
 *  @param va <#va description#>
 *  @param vn <#vn description#>
 */
void
VECTOR4D_Normalize(VECTOR4D_PTR va,VECTOR4D_PTR vn) {
    VECTOR4D_ZERO(vn);
    // compute length
    float length = sqrt(va->x*va->x + va->y*va->y + va->z*va->z);
    
    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;
    
    float length_inv = 1.0/length;
    
    // compute normalized version of vector
    vn->x = va->x*length_inv;
    vn->y = va->y*length_inv;
    vn->z = va->z*length_inv;
    vn->w = 1;
}
/**
 these are the 4D version of the vector functions, they assume that the vectors
  are 3D with a w, so w is left ;out of all the operations
 */
void
VECTOR4D_Build(VECTOR4D_PTR init ,VECTOR4D_PTR term , VECTOR4D_PTR result) {
    // build a 4d vector
    result->x = term->x - init->x;
    result->y = term->y - init->y;
    result->z = term->z - init->z;
    result->w = 1;
}
/**
 *  this function returns the cosine of the angle between,two vectors. Note, 
 *  we could compute the actual angle ;many many times, in further calcs we
 *  will want ultimately compute cos of the angle, so why not just leave it!
 *
 *  @param va <#va description#>
 *  @param vb <#vb description#>
 *
 *  @return <#return value description#>
 */
float
VECTOR4D_CosTh(VECTOR4D_PTR va, VECTOR4D_PTR vb) {
    return(VECTOR4D_Dot(va,vb)/(VECTOR4D_Length(va)*VECTOR4D_Length(vb)));
}

void
VECTOR4D_Print(VECTOR4D_PTR va, char *name) {
    Write_Error("\n%s[",name);
    for (int index=0; index<4; index++)
        Write_Error("%f, ",va->M[index]);
    Write_Error("]");

}

/******************************************************************************/

// 2x2 matrix functions (note there others in T3DLib1.cpp|h)
/**
 *  this function fills a 2X2 matrix with the sent data in row major form
 *
 *  @param ma  <#ma description#>
 *  @param m00 <#m00 description#>
 *  @param m01 <#m01 description#>
 *  @param m10 <#m10 description#>
 *  @param m11 <#m11 description#>
 */
void
Mat_Init_2X2(MATRIX2X2_PTR ma,float m00,float m01,float m10,float m11) {
    ma->M00 = m00; ma->M01 = m01;ma->M10=m10; ma->M11=m11;
}

void
Print_Mat_2X2(MATRIX2X2_PTR ma, char *name) {
    // prints out a 3x3 matrix
    Write_Error("\n%s=\n",name);
    
    for (int r=0; r < 2; r++, Write_Error("\n"))
        for (int c=0; c < 2; c++)
            Write_Error("%f ",ma->M[r][c]);
}

/** 用于计算矩阵m的行列式,并将结果返回到堆栈
 *  computers the determinate of a 2X2 matrix
 *
 *  @param m <#m description#>
 *
 *  @return <#return value description#>
 */
float
Mat_Det_2X2(MATRIX2X2_PTR m) {
    return (m->M00 * m->M11 - m->M01 * m->M10);
}
/**
 *  this functions adds two 2X2 matrices together and stores the result in msum
 *
 *  @param ma   <#ma description#>
 *  @param mb   <#mb description#>
 *  @param msum <#msum description#>
 */
void
Mat_Add_2X2(MATRIX2X2_PTR ma,MATRIX2X2_PTR mb,MATRIX2X2_PTR msum) {
    msum->M00 = ma->M00 + mb->M00;
    msum->M01 = ma->M01 + mb->M01;
    msum->M10 = ma->M10+mb->M10;
    msum->M11 = ma->M11+mb->M11;
}
/**
 *  this function multiplies two 2X2 matrices together and stores the result in mprod
 *
 *  @param ma    <#ma description#>
 *  @param mb    <#mb description#>
 *  @param mprod <#mprod description#>
 */
void
Mat_Mul_2X2(MATRIX2X2_PTR ma,MATRIX2X2_PTR mb,MATRIX2X2_PTR mprod) {
    mprod->M00 = ma->M00*mb->M00 + ma->M01*mb->M10;
    mprod->M01 = ma->M00*mb->M01 + ma->M01*mb->M11;
    
    mprod->M10 = ma->M10*mb->M00 + ma->M11*mb->M10;
    mprod->M11 = ma->M10*mb->M01 + ma->M11*mb->M11;
}
/** 用于计算矩阵m的逆矩阵(如果有的情况),并将结果存储在mi中,若逆矩阵存在则返回1,否则返回0且mi为定义
 *  this function computers the inverse of a 2x2 matrix and stores the result in mi
 *
 *  @param ma <#ma description#>
 *  @param mi <#mi description#>
 *
 *  @return <#return value description#>
 */
int
Mat_Inverse_2X2(MATRIX2X2_PTR ma,MATRIX2X2_PTR mi) {
    // compute determinate
    float det = (ma->M00*ma->M11 - ma->M01*ma->M10);
    
    // if determinate is 0 then inverse doesn't exist
    if (fabs(det) < EPSILON_E5) return 0;
    
    float det_inv = 1.0/det;
    // fill in inverse by formula
    mi->M00 =  ma->M11*det_inv;
    mi->M01 = -ma->M01*det_inv;
    mi->M10 = -ma->M10*det_inv;
    mi->M11 =  ma->M00*det_inv;
    // return sucess
    return 1;
}
/**
 *  solves the system AX=B and computes X=A(-1)*B
 *  by using cramers rule and determinates
 *  @param A <#A description#>
 *  @param X <#X description#>
 *  @param B <#B description#>
 *
 *  @return <#return value description#>
 */
int
Solve_2X2_System(MATRIX2X2_PTR A,MATRIX1X2_PTR X,MATRIX1X2_PTR B) {
    // step 1: compute determinate of A
    float det_A = Mat_Det_2X2(A);
    
    // test if det(a) is zero, if so then there is no solution
    if (fabs(det_A) < EPSILON_E5)
        return(0);
    
    // step 2: create x,y numerator matrices by taking A and
    // replacing each column of it with B(transpose) and solve
    MATRIX2X2 work_mat; // working matrix
    
    // solve for x /////////////////
    
    // copy A into working matrix
    MAT_COPY_2X2(A, &work_mat);
    
    // swap out column 0 (x column)
    MAT_COLUMN_SWAP_2X2(&work_mat, 0, B);
    
    // compute determinate of A with B swapped into x column
    float det_ABx = Mat_Det_2X2(&work_mat);
    
    // now solve for X00
    X->M00 = det_ABx/det_A;
    
    // solve for y /////////////////
    
    // copy A into working matrix
    MAT_COPY_2X2(A, &work_mat);
    
    // swap out column 1 (y column)
    MAT_COLUMN_SWAP_2X2(&work_mat, 1, B);
    
    // compute determinate of A with B swapped into y column
    float det_ABy = Mat_Det_2X2(&work_mat);
    
    // now solve for X01
    X->M01 = det_ABy/det_A;
    
    // return success
    return 1;
}
/******************************************************************************/


//3X3 matrix functions
/** 将一个1X2矩阵(实质就是一个2D点)与一个3X3矩阵(表示旋转或平移)相乘.这种运算在数学上时未定义的,
 * 因为他们的内纬不相同,然而,如果假设1X2实际是1X3矩阵,其最后一个元素为1,将可以执行这种乘法,
 *  这个函数对2D点进行变换时 非常有用。
 *  waring:这个函数速度很慢,可以消除循环来加速.但现在更清晰,必要时可以优化，对代码优化之前应该
 * 先对算法进行优化。有了最优的算法和代码后,可以有内联汇编
 *  this function mutiplies a 1X2 matrix against a 3X2 matrix - ma*mb and stores the
 *  result using a dummy element for the 3rd element of the 1X2 to make the matrix
 *  multiply vaild . 1X3 X 3X2
 *  @param ma    <#ma description#>
 *  @param mb    <#mb description#>
 *  @param mprod <#mprod description#>
 *
 *  @return <#return value description#>
 */
int
MAT_Mul_1X2_3X2(MATRIX1X2_PTR ma, MATRIX3X2_PTR mb,MATRIX1X2_PTR mprod) {
    for (int col = 0; col < 2; col++) {
        //computer dot product from row of ma  and  column of mb
        float sum =0;//used to hold result
        int index;
        for (index = 0; index < 2; index++) {
            // add in next product pair
            sum+=(ma->M[index]*mb->M[index][col]);
        }
        //add in last element * 1
        sum += mb->M[index][col];
        
        //insert resulting col  element
        mprod->M[col] = sum;
    }
    return 1;
}
/** 将一个1x3矩阵(实质是一个行向量)与一个3x3矩阵相乘.等价于Mat_Mul_VECTOR3D_3X3(),只是类型不同
 *  this function multiplies a 1X3 matrix against a 3X3 matrix - ma*mb and stores
 *   the result
 *  @param ma    <#ma description#>
 *  @param mb    <#mb description#>
 *  @param mprod <#mprod description#>
 *
 *  @return <#return value description#>
 */
int
Mat_Mul_1X3_3X3(MATRIX1X3_PTR ma, MATRIX3X3_PTR mb,MATRIX1X3_PTR mprod) {
    for (int col = 0; col < 3; col++) {
        float sum = 0;// used to hold result
        int index;
        for (index = 0; index < 3 ; index++) {
            sum += (ma->M[index]*mb->M[index][col]);
        }
        //insert resulting col element
        mprod->M[col] =sum;
    }
    return 1;
}
/** 将两个矩阵相乘(ma * mb),并将结果存储到mprod中
 *  this function multiplies two matrices together and stores the result
 *
 *  @param ma    <#ma description#>
 *  @param mb    <#mb description#>
 *  @param mprod <#mprod description#>
 *
 *  @return <#return value description#>
 */
int
Mat_Mul_3X3(MATRIX3X3_PTR ma,MATRIX3X3_PTR mb,MATRIX3X3_PTR mprod) {
    for (int row = 0; row < 3; row++) {
        for (int col = 0 ; col < 3; col++) {
            float sum = 0; //used to hold result
            for (int index = 0; index < 3; index++) {
                sum += (ma->M[row][index] *mb->M[index][col]);
            }
            //insert resulting row col element
            mprod->M[row][col] =sum;
        }
    }
    return 1;
}

/** 使用传入的浮点值以先行后列的次序 初始化矩阵ma
 *  this function fills a  3X2 matrix with the sent data in row major form
 *
 *  @param ma  <#ma description#>
 *  @param m00 <#m00 description#>
 *  @param m01 <#m01 description#>
 *  @param m10 <#m10 description#>
 *  @param m11 <#m11 description#>
 *  @param m20 <#m20 description#>
 *  @param m21 <#m21 description#>
 *
 *  @return <#return value description#>
 */
inline int
Mat_Init_3X2(MATRIX3X2_PTR ma,float m00, float m01,float m10, float m11,
             float m20, float m21) {
    ma->M[0][0] = m00;ma->M[0][1] = m01;ma->M[1][0]=m01;ma->M[1][1]=m11;
    
    return 1;
}

/** 将两个矩阵相加(ma＋mb),并将结果存储到msum中
 *  this function adds two 3X3 matrices together and stores the result.
 *
 *  @param ma   <#ma description#>
 *  @param mb   <#mb description#>
 *  @param msum <#msum description#>
 */
void
Mat_Add_3X3(MATRIX3X3_PTR ma, MATRIX3X3_PTR mb, MATRIX3X3_PTR msum) {
    for (int row = 0; row < 3; row++) {
        for (int col = 0 ; col < 3; col++) {
            msum->M[row][col] = ma->M[row][col] + mb->M[row][col];
        }
    }
}
/** 将1x3 的行向量(3d点) 与 3x3矩阵相乘,并将结果存储到1X3的行向量vprod (实质 执行点或向量 与矩阵的乘法)
 *  this function multiplies a VECTOR3D against a 3X3 matrix - ma*mb and
 *  stores the result in vprod
 *  @param va    <#va description#>
 *  @param mb    <#mb description#>
 *  @param vprod <#vprod description#>
 */
void
Mat_Mul_VECTOR3D_3X3(VECTOR3D_PTR  va, MATRIX3X3_PTR mb,VECTOR3D_PTR  vprod) {
    for (int col = 0 ; col < 3; col++) {
        float sum = 0;
        for (int row = 0; row < 3; row++) {
            sum+=(va->M[row]*mb->M[row][col]);
        }
        vprod->M[col] = sum;
    }
}
/** 使用传入的浮点值以先行后列的次序 初始化矩阵ma
 *  this function fills a 3x3 matrix with sent data in row major form
 *
 *  @param ma  <#ma description#>
 *  @param m00 <#m00 description#>
 *  @param m01 <#m01 description#>
 *  @param m02 <#m02 description#>
 *  @param m10 <#m10 description#>
 *  @param m11 <#m11 description#>
 *  @param m12 <#m12 description#>
 *  @param m20 <#m20 description#>
 *  @param m21 <#m21 description#>
 *  @param m22 <#m22 description#>
 */
void
Mat_Init_3X3(MATRIX3X3_PTR ma,float m00, float m01, float m02,float m10, float m11, float m12,
             float m20, float m21, float m22) {
    ma->M00 = m00; ma->M01 = m01; ma->M02 = m02;
    ma->M10 = m10; ma->M11 = m11; ma->M12 = m12;
    ma->M20 = m20; ma->M21 = m21; ma->M22 = m22;
}


/** 用于计算矩阵m的逆矩阵(如果有的情况),并将结果存储在mi中,若逆矩阵存在则返回1,否则返回0且mi为定义
 *  this function computers the inverse of a 3x3
  *
 *  @param m  <#m description#>
 *  @param mi <#mi description#>
 *
 *  @return <#return value description#>
 */
int
Mat_Inverse_3X3(MATRIX3X3_PTR m, MATRIX3X3_PTR mi) {
    //first compute the determinate to see if there is an inverse
    float det = m->M00*(m->M11*m->M22 - m->M21*m->M12) -
    m->M01*(m->M10*m->M22 - m->M20*m->M12) +
    m->M02*(m->M10*m->M21 - m->M20*m->M11);
    
    if (fabs(det) < EPSILON_E5)
        return(0);
    
    // compute inverse to save divides
    float det_inv = 1.0/det;
    
    // compute inverse using m-1 = adjoint(m)/det(m)
    mi->M00 =  det_inv*(m->M11*m->M22 - m->M21*m->M12);
    mi->M10 = -det_inv*(m->M10*m->M22 - m->M20*m->M12);
    mi->M20 =  det_inv*(m->M10*m->M21 - m->M20*m->M11);
    
    mi->M01 = -det_inv*(m->M01*m->M22 - m->M21*m->M02);
    mi->M11 =  det_inv*(m->M00*m->M22 - m->M20*m->M02);
    mi->M21 = -det_inv*(m->M00*m->M21 - m->M20*m->M01);
    
    mi->M02 =  det_inv*(m->M01*m->M12 - m->M11*m->M02);
    mi->M12 = -det_inv*(m->M00*m->M12 - m->M10*m->M02);
    mi->M22 =  det_inv*(m->M00*m->M11 - m->M10*m->M01);
    
    return 1;
}
/** 计算矩阵m的行列式,并将结果返回到堆栈
 *  computes the determinate of a 3x3 matrix using co-factor expansion
 *
 *  @param m <#m description#>
 *
 *  @return <#return value description#>
 */
float
Mat_Det_3X3(MATRIX3X3_PTR m) {
    return(m->M00*(m->M11*m->M22 - m->M21*m->M12) -
           m->M01*(m->M10*m->M22 - m->M20*m->M12) +
           m->M02*(m->M10*m->M21 - m->M20*m->M11) );

}
/**
 *  solves the system AX=B and computes X=A(-1)*B by using cramers rule and determinates
 *
 *  @param A <#A description#>
 *  @param X <#X description#>
 *  @param B <#B description#>
 *
 *  @return <#return value description#>
 */
int
Solve_3X3_System(MATRIX3X3_PTR A, MATRIX1X3_PTR X, MATRIX1X3_PTR B) {
    // step 1: compute determinate of A
    float det_A = Mat_Det_3X3(A);
    
    // test if det(a) is zero, if so then there is no solution
    if (fabs(det_A) < EPSILON_E5)
        return(0);
    
    // step 2: create x,y,z numerator matrices by taking A and
    // replacing each column of it with B(transpose) and solve
    MATRIX3X3 work_mat; // working matrix
    
    // solve for x /////////////////
    
    // copy A into working matrix
    MAT_COPY_3X3(A, &work_mat);
    
    // swap out column 0 (x column)
    MAT_COLUMN_SWAP_3X3(&work_mat, 0, B);
    
    // compute determinate of A with B swapped into x column
    float det_ABx = Mat_Det_3X3(&work_mat);
    
    // now solve for X00
    X->M00 = det_ABx/det_A;
    
    // solve for y /////////////////
    
    // copy A into working matrix
    MAT_COPY_3X3(A, &work_mat);
    
    // swap out column 1 (y column)
    MAT_COLUMN_SWAP_3X3(&work_mat, 1, B);
    
    // compute determinate of A with B swapped into y column
    float det_ABy = Mat_Det_3X3(&work_mat);
    
    // now solve for X01
    X->M01 = det_ABy/det_A;
    
    // solve for z /////////////////
    
    // copy A into working matrix
    MAT_COPY_3X3(A, &work_mat);
    
    // swap out column 2 (z column)
    MAT_COLUMN_SWAP_3X3(&work_mat, 2, B);
    
    // compute determinate of A with B swapped into z column
    float det_ABz = Mat_Det_3X3(&work_mat);
    
    // now solve for X02
    X->M02 = det_ABz/det_A;
    return 1;
}
/******************************************************************************/






/******************************************************************************/







/******************************************************************************/





/******************************************************************************/






