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
//4x4 matrix functions
/**
 *  this  function adds two 4x4 matrices together and stores the result
 *
 *  @param ma   <#ma description#>
 *  @param mb   <#mb description#>
 *  @param msum <#msum description#>
 */
void
Mat_Add_4X4(MATRIX4X4_PTR ma,MATRIX4X4_PTR mb,MATRIX4X4_PTR msum) {
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            msum->M[row][col] = ma->M[row][col] + mb->M[row][col];            
        }
    }
}
/**
 *  this function multiplies two 4x4 matrices together and stores the result in mprod
 *  note later we will take advantage of the fact that we know that w=1 always, 
 *  and that the last column of a 4x4 is always 0
 *  @param ma    <#ma description#>
 *  @param mb    <#mb description#>
 *  @param mprod <#mprod description#>
 */
void
Mat_Mul_4X4(MATRIX4X4_PTR ma,MATRIX4X4_PTR mb,MATRIX4X4_PTR mprod) {
    for (int row =0 ; row < 4; row++) {
        for (int col = 0; col<4; col++) {
            float sum =0 ;
            for (int index = 0; index < 4; index++) {
                sum+=(ma->M[row][index]*mb->M[index][col]);
            }
            mprod->M[row][col] = sum;
        }
    }
}
/**
 *  this function multiplies a 1x4 matrix against a 4x4 matrix - ma*mb and stores the result
 *  no tricks or assumptions here, just a straight multiply
 *  @param ma    <#ma description#>
 *  @param mb    <#mb description#>
 *  @param mprod <#mprod description#>
 */
void
Mat_Mul_1X4_4X4(MATRIX1X4_PTR ma,MATRIX4X4_PTR mb,MATRIX1X4_PTR mprod) {
    for (int col=0; col<4; col++) {
        float sum = 0; //used to hold result
        for (int row = 0 ; row< 4; row++) {
            sum+=(ma->M[row] * mb->M[row][col]);
        }
        mprod->M[col] = sum;
    }
}
/**
 *  this function multiplies a VECTOR3D against a 4x4 matrix - ma*mb and stores the result in mprod
 *  the function assumes that the vector refers to a 4D homogenous vector, thus the function assumes that
 *  w=1 to carry out the multiply, also the function does not carry out the last column multiply since
 *  we are assuming w=1, there is no point
 *  @param va    <#va description#>
 *  @param mb    <#mb description#>
 *  @param mprod <#mprod description#>
 */
void
Mat_Mul_VECTOR3D_4X4(VECTOR3D_PTR  va,MATRIX4X4_PTR mb,VECTOR3D_PTR  vprod) {
    for (int col=0; col < 3; col++) {
        // compute dot product from row of ma and column of mb
        float sum = 0; // used to hold result
        int row ;
        for (row =0; row<3; row++) {
            // add in next product pair
            sum+=(va->M[row]*mb->M[row][col]);
        }
        sum += mb->M[row][col];
        // insert resulting col element
        vprod->M[col] = sum;
    }
}
/**
 *   this function multiplies a VECTOR3D against a
 *  4x3 matrix - ma*mb and stores the result in mprod
 *  the function assumes that the vector refers to a
 *  4D homogenous vector, thus the function assumes that
 *  w=1 to carry out the multiply, also the function
 *  does not carry out the last column multiply since
 *  we are assuming w=1, there is no point
 *
 *  @param va    <#va description#>
 *  @param mb    <#mb description#>
 *  @param vprod <#vprod description#>
 */
void
Mat_Mul_VECTOR3D_4X3(VECTOR3D_PTR  va,MATRIX4X3_PTR mb,VECTOR3D_PTR  vprod) {
    for (int col=0; col < 3; col++){
        float sum = 0; // used to hold result
        int row;
        for (row=0; row<3; row++){
            sum+=(va->M[row]*mb->M[row][col]);
        }
        sum+=mb->M[row][col];
        // insert resulting col element
        vprod->M[col] = sum;
    }
}
/**
 *  this function multiplies a VECTOR4D against a
 *  4x4 matrix - ma*mb and stores the result in mprod
 *  the function makes no assumptions
 *
 *  @param va    <#va description#>
 *  @param mb    <#mb description#>
 *  @param vprod <#vprod description#>
 */
void
Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR  va,MATRIX4X4_PTR mb,VECTOR4D_PTR  vprod) {
    for (int col=0; col < 4; col++){
        float sum = 0; // used to hold result
        int row;
        for (row=0; row<4; row++){
            sum+=(va->M[row]*mb->M[row][col]);
        }
        vprod->M[col] = sum;
    }
}
/**
 *  this function multiplies a VECTOR4D against a
 *  4x3 matrix - ma*mb and stores the result in mprod
 *  the function assumes that the last column of
 *  mb is [0 0 0 1]t , thus w just gets replicated
 *  from the vector [x y z w]
 *
 *  @param va    <#va description#>
 *  @param mb    <#mb description#>
 *  @param vprod <#vprod description#>
 */
void
Mat_Mul_VECTOR4D_4X3(VECTOR4D_PTR  va,MATRIX4X4_PTR mb,VECTOR4D_PTR  vprod) {
    for (int col=0; col < 3; col++) {
        float sum = 0; // used to hold result
        for (int row=0; row<4; row++){
            sum+=(va->M[row]*mb->M[row][col]);
        }
        vprod->M[col] = sum;
    }
    vprod->M[3] = va->M[3];
}
/**
 *  this function fills a 4x4 matrix with the sent data in row major form
 */
void
Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33) {
    
    ma->M00 = m00; ma->M01 = m01; ma->M02 = m02; ma->M03 = m03;
    ma->M10 = m10; ma->M11 = m11; ma->M12 = m12; ma->M13 = m13;
    ma->M20 = m20; ma->M21 = m21; ma->M22 = m22; ma->M23 = m23;
    ma->M30 = m30; ma->M31 = m31; ma->M32 = m32; ma->M33 = m33;
}
/**
 *  computes the inverse of a 4x4, assumes that the last column is [0 0 0 1]t
 *
 *  @param m  <#m description#>
 *  @param mi <#mi description#>
 *
 *  @return <#return value description#>
 */
int
Mat_Inverse_4X4(MATRIX4X4_PTR m, MATRIX4X4_PTR mi) {
    float det =  ( m->M00 * ( m->M11 * m->M22 - m->M12 * m->M21 ) -
                  m->M01 * ( m->M10 * m->M22 - m->M12 * m->M20 ) +
                  m->M02 * ( m->M10 * m->M21 - m->M11 * m->M20 ) );
    // test determinate to see if it's 0
    if (fabs(det) < EPSILON_E5) return 0;
    
    float det_inv  = 1.0f / det;
    mi->M00 =  det_inv * ( m->M11 * m->M22 - m->M12 * m->M21 );
    mi->M01 = -det_inv * ( m->M01 * m->M22 - m->M02 * m->M21 );
    mi->M02 =  det_inv * ( m->M01 * m->M12 - m->M02 * m->M11 );
    mi->M03 = 0.0f; // always 0
    
    mi->M10 = -det_inv * ( m->M10 * m->M22 - m->M12 * m->M20 );
    mi->M11 =  det_inv * ( m->M00 * m->M22 - m->M02 * m->M20 );
    mi->M12 = -det_inv * ( m->M00 * m->M12 - m->M02 * m->M10 );
    mi->M13 = 0.0f; // always 0
    
    mi->M20 =  det_inv * ( m->M10 * m->M21 - m->M11 * m->M20 );
    mi->M21 = -det_inv * ( m->M00 * m->M21 - m->M01 * m->M20 );
    mi->M22 =  det_inv * ( m->M00 * m->M11 - m->M01 * m->M10 );
    mi->M23 = 0.0f; // always 0
    
    mi->M30 = -( m->M30 * mi->M00 + m->M31 * mi->M10 + m->M32 * mi->M20 );
    mi->M31 = -( m->M30 * mi->M01 + m->M31 * mi->M11 + m->M32 * mi->M21 );
    mi->M32 = -( m->M30 * mi->M02 + m->M31 * mi->M12 + m->M32 * mi->M22 );
    mi->M33 = 1.0f; // always 0
    return 1;
}
/******************************************************************************/
//四元数函数
/**
 *  this function adds two quaternions
 *
 *  @param q1   <#q1 description#>
 *  @param q2   <#q2 description#>
 *  @param qsum <#qsum description#>
 */
void
QUAT_Add(QUAT_PTR q1,QUAT_PTR q2,QUAT_PTR qsum){
    qsum->x = q1->x + q2->x;
    qsum->y = q1->y + q2->y;
    qsum->z = q1->z + q2->z;
    qsum->w = q1->w + q2->w;
}
/**
 *  this function subtracts two quternions; (q1 - q2)
 *
 *  @param q1    <#q1 description#>
 *  @param q2    <#q2 description#>
 *  @param qdiff <#qdiff description#>
 */
void
QUAT_Sub(QUAT_PTR q1,QUAT_PTR q2,QUAT_PTR qdiff) {
    qdiff->x = q1->x - q2->x;
    qdiff->y = q1->y - q2->y;
    qdiff->z = q1->z - q2->z;
    qdiff->w = q1->w - q2->w;
}
/**
 *  this function computers the conjugate of a quterninons
 *
 *  @param q     <#q description#>
 *  @param qconj <#qconj description#>
 */
void
QUAT_Conjugate(QUAT_PTR q,QUAT_PTR qconj) {
    qconj->x = -q->x;
    qconj->y = -q->y;
    qconj->z = -q->z;
    qconj->w = q->w;
}
/**
 *  this function scales a quaternion and return it
 *
 *  @param q     <#q description#>
 *  @param scale <#scale description#>
 *  @param qs    <#qs description#>
 */
void
QUAT_Scale(QUAT_PTR q,float scale,QUAT_PTR qs) {
    qs->x = scale * q->x;
    qs->y = scale*q->y;
    qs->z = scale*q->z;
    qs->w = scale*q->w;
}
/**
 *  this function scales a quaternion in place
 *
 *  @param q     <#q description#>
 *  @param scale <#scale description#>
 */
void
QUAT_Scale(QUAT_PTR q,float scale) {
    q->x *= scale;
    q->y *= scale;
    q->z *= scale;
    q->w *= scale;

}
/**
 *  returns the length or norm of a quaterion
 *
 *  @param q <#q description#>
 *
 *  @return <#return value description#>
 */
float
QUAT_Norm(QUAT_PTR q) {
    return(sqrtf(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z));
}
/**
 *  returns the length or norm of a quaternion squared
 *
 *  @param q <#q description#>
 *
 *  @return <#return value description#>
 */
float
QUAT_Norm2(QUAT_PTR q) {
    return(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
}
/**
 *  this functions normalizes the sent quaternion and returns it
 *
 *  @param q  <#q description#>
 *  @param qn <#qn description#>
 */
void
QUAT_Normalize(QUAT_PTR q,QUAT_PTR qn) {
    // compute 1/length
    float qlength_inv = 1.0/(sqrtf(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z));
    
    // now normalize
    qn->w=q->w*qlength_inv;
    qn->x=q->x*qlength_inv;
    qn->y=q->y*qlength_inv;
    qn->z=q->z*qlength_inv;
}
/**
 *  this functions normalizes the sent quaternion in place
 *
 *  @param q <#q description#>
 */
void
QUAT_Normalize(QUAT_PTR q) {
    // compute length
    float qlength_inv = 1.0/(sqrtf(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z));
    // now normalize
    q->w*=qlength_inv;
    q->x*=qlength_inv;
    q->y*=qlength_inv;
    q->z*=qlength_inv;

}
/**
 *  this function computers the inverse of a unit quaternion and returns the result
 *
 *  @param q  <#q description#>
 *  @param qi <#qi description#>
 */
void
QUAT_Unit_Inverse(QUAT_PTR q, QUAT_PTR qi) {
    // the inverse of a unit quaternion is the conjugate :)
    qi->w =  q->w;
    qi->x = -q->x;
    qi->y = -q->y;
    qi->z = -q->z;
}
/**
 *  this function computes the inverse of a unit quaternion in place
 *
 *  @param q <#q description#>
 */
void
QUAT_Unit_Inverse(QUAT_PTR q) {
    // the inverse of a unit quaternion is the conjugate :)
    q->x = -q->x;
    q->y = -q->y;
    q->z = -q->z;
}
/**
 *  this function computers the inverse of general quaternion and returns result
 *  in general, q-1 = *q/|q|2
 *  @param q  <#q description#>
 *  @param qi <#qi description#>
 */
void
QUAT_Inverse(QUAT_PTR q, QUAT_PTR qi) {
    // compute norm squared
    float norm2_inv = 1.0/(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
    
    // and plug in
    qi->w =  q->w*norm2_inv;
    qi->x = -q->x*norm2_inv;
    qi->y = -q->y*norm2_inv;
    qi->z = -q->z*norm2_inv;
}
/**
 *  this function computers the inverse of a general quaternion in place
 *
 *  @param q <#q description#>
 */
void
QUAT_Inverse(QUAT_PTR q) {
    // in general, q-1 = *q/|q|2
    // compute norm squared
    float norm2_inv = 1.0/(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
    
    // and plug in
    q->w =  q->w*norm2_inv;
    q->x = -q->x*norm2_inv;
    q->y = -q->y*norm2_inv;
    q->z = -q->z*norm2_inv;

}
/**
 *  this function mutiplies two quaternions
 *
 *  @param q1    <#q1 description#>
 *  @param q2    <#q2 description#>
 *  @param qprod <#qprod description#>
 */
void
QUAT_Mul(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qprod) {
    
    // this is the brute force method
    //qprod->w = q1->w*q2->w - q1->x*q2->x - q1->y*q2->y - q1->z*q2->z;
    //qprod->x = q1->w*q2->x + q1->x*q2->w + q1->y*q2->z - q1->z*q2->y;
    //qprod->y = q1->w*q2->y - q1->x*q2->z + q1->y*q2->w - q1->z*q2->x;
    //qprod->z = q1->w*q2->z + q1->x*q2->y - q1->y*q2->x + q1->z*q2->w;
    
    // this method was arrived at basically by trying to factor the above
    // expression to reduce the # of multiplies
    
    float prd_0 = (q1->z - q1->y) * (q2->y - q2->z);
    float prd_1 = (q1->w + q1->x) * (q2->w + q2->x);
    float prd_2 = (q1->w - q1->x) * (q2->y + q2->z);
    float prd_3 = (q1->y + q1->z) * (q2->w - q2->x);
    float prd_4 = (q1->z - q1->x) * (q2->x - q2->y);
    float prd_5 = (q1->z + q1->x) * (q2->x + q2->y);
    float prd_6 = (q1->w + q1->y) * (q2->w - q2->z);
    float prd_7 = (q1->w - q1->y) * (q2->w + q2->z);
    
    float prd_8 = prd_5 + prd_6 + prd_7;
    float prd_9 = 0.5 * (prd_4 + prd_8);
    
    // and finally build up the result with the temporary products
    
    qprod->w = prd_0 + prd_9 - prd_5;
    qprod->x = prd_1 + prd_9 - prd_8;
    qprod->y = prd_2 + prd_9 - prd_7;
    qprod->z = prd_3 + prd_9 - prd_6;

}
/**
 *  this function computers q1* q2*q3 in that order and returns the results in qprod
 *
 *  @param q1    <#q1 description#>
 *  @param q2    <#q2 description#>
 *  @param q3    <#q3 description#>
 *  @param qprod <#qprod description#>
 */
void
QUAT_Triple_Product(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR q3,
                         QUAT_PTR qprod) {
    QUAT qtmp;
    QUAT_Mul(q1,q2,&qtmp);
    QUAT_Mul(&qtmp, q3, qprod);
}
/**
 *  initializes a quaternion based on a 3d direction vector and angle
 *  note the direction vector must be a unit vector and the angle is in rads
 *  @param q     <#q description#>
 *  @param v     <#v description#>
 *  @param theta <#theta description#>
 */
void
VECTOR3D_Theta_To_QUAT(QUAT_PTR q, VECTOR3D_PTR v, float theta) {
    float theta_div_2 = (0.5)*theta; // compute theta/2
    // pre-compute to save time
    float sinf_theta = sinf(theta_div_2);
    
    q->x = sinf_theta * v->x;
    q->y = sinf_theta * v->y;
    q->z = sinf_theta * v->z;
    q->w = cosf(theta_div_2);
}
/**
 *  initializes a quaternion based on a 4d direction vector and angle
 *  note the direction vector must be a unit vector and the angle is in rads
 *  @param q     <#q description#>
 *  @param v     <#v description#>
 *  @param theta <#theta description#>
 */
void
VECTOR4D_Theta_To_QUAT(QUAT_PTR q, VECTOR4D_PTR v, float theta) {
    float theta_div_2 = (0.5)*theta; // compute theta/2
    // pre-compute to save time
    float sinf_theta = sinf(theta_div_2);
    
    q->x = sinf_theta * v->x;
    q->y = sinf_theta * v->y;
    q->z = sinf_theta * v->z;
    q->w = cosf(theta_div_2);
}
/**
 *  this function intializes a quaternion based on the zyx multiplication order
 *   of the angles that are parallel to the zyx axis respectively.
 *   note there are 11 other possibilities
 *
 *  @param q       <#q description#>
 *  @param theta_z <#theta_z description#>
 *  @param theta_y <#theta_y description#>
 *  @param theta_x <#theta_x description#>
 */
void
EulerZYX_To_QUAT(QUAT_PTR q, float theta_z, float theta_y, float theta_x) {
    // precompute values
    float cos_z_2 = 0.5*cosf(theta_z);
    float cos_y_2 = 0.5*cosf(theta_y);
    float cos_x_2 = 0.5*cosf(theta_x);
    
    float sin_z_2 = 0.5*sinf(theta_z);
    float sin_y_2 = 0.5*sinf(theta_y);
    float sin_x_2 = 0.5*sinf(theta_x);
    
    // and now compute quaternion
    q->w = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
    q->x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
    q->y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
    q->z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;
}
/**
 *  this function converts a unit quaternion into a unit direction
 *  vector and rotation angle about that vector
 *  @param q     <#q description#>
 *  @param v     <#v description#>
 *  @param theta <#theta description#>
 */
void
QUAT_To_VECTOR3D_Theta(QUAT_PTR q, VECTOR3D_PTR v, float *theta) {
    // extract theta
    *theta = acosf(q->w);
    // pre-compute to save time
    float sinf_theta_inv = 1.0/sinf(*theta);
    // now the vector
    v->x    = q->x*sinf_theta_inv;
    v->y    = q->y*sinf_theta_inv;
    v->z    = q->z*sinf_theta_inv;
    // multiply by 2
    *theta*=2;
}
/******************************************************************************/

//2D参数化直线函数
/** 根据指定的点计算它们之间的向量,并初始化一条2D参数化直线,note:该向量是在函数内部生成
 *  this initializes a parametric 2d line
 *  note that the function computes v=p_pend - p_init, thus when t=0 the line p=p0+v*t = p0
 *  and whne t=1, p=p1, this way the segement is traced out from p0 to p1 via 0<= t <= 1
 */
void
Init_Parm_Line2D(POINT2D_PTR p_init,POINT2D_PTR p_term,PARMLINE2D_PTR p) {
    // start point
    VECTOR2D_INIT(&(p->p0), p_init);
    // end point
    VECTOR2D_INIT(&(p->p1), p_term);
    // now compute direction vector from p0->p1
    VECTOR2D_Build(p_init, p_term, &(p->v));
}

/**  计算2D参数化直线在参数t处的值,并将其返回到指定的点中 ;当t＝0,结果为起点p1,当t=1,结果为终点p2,也就是说当t从0变化到1,将从p1移动到p2.
 *  this function computers the value of the sent parametric line at the value of t
 *
 *  @param p  <#p description#>
 *  @param t  <#t description#>
 *  @param pt <#pt description#>
 */
void
Compute_Parm_Line2D(PARMLINE2D_PTR p, float t, POINT2D_PTR pt) {
    pt->x = p->p0.x + p->v.x*t;
    pt->y = p->p0.y + p->v.y*t;
}
/** waring:该函数不检测两条直线是否共线,因为这将极大地降低速度,同时有太多的情况(部分重叠,包含,只有一点重叠等情况)。
 *  这个函数计算两条参数化线段的交点,并将t1和t2分别设置交点在p1和p2上对应的t值
 *  然而t值可能不在范围[0,1] 内,这意味着线段本身没有相交,但它们对应的直线是相交的
 *  函数返回0表示没有相交,返回1表示交点在线段上 ;2代表相交,但交点不在线段上,3代表两条线段位于同一直线上
 */
int
Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2,
                           float *t1, float *t2) {
    
    // step 1: 检测它们是否平行
    // 如果一个方向向量是另一个向量与一个标量的乘积,则说明两条线段平行
    // 除非它们重叠,否则不可能相交
    float det_p1p2 = (p1->v.x*p2->v.y - p1->v.y*p2->v.x);
    if (fabs(det_p1p2) <= EPSILON_E5){
        //这表明两条线段要么根本不相交,要么位于同一条直线上
        //在后一种情况下,可能有一个或多个交点
        //现在暂时假设它们不相交,以后需要时再重编写该函数,以考虑重叠的情况
        return(PARM_LINE_NO_INTERSECT);
    }
    // step 2: 计算t1和t2的值;我们有两条以下述方式表示的线段
    // p    = p0    +  v*t, specifically
    // p1   = p10   + v1*t1
    // p1.x = p10.x + v1.x*t1
    // p1.y = p10.y + v1.y*t1
    
    // p2 = p20 + v2*t2
    // p2   = p20   + v2*t2
    // p2.x = p20.x + v2.x*t2
    // p2.y = p20.y + v2.y*t2
    // solve the system when x1 = x2 and y1 = y2
    // 计算交点
    *t1 = (p2->v.x*(p1->p0.y - p2->p0.y) - p2->v.y*(p1->p0.x - p2->p0.x)) /det_p1p2;
    *t2 = (p1->v.x*(p1->p0.y - p2->p0.y) - p1->v.y*(p1->p0.x - p2->p0.x)) /det_p1p2;
    // 检测交点是否在线段上
    if ((*t1>=0) && (*t1<=1) && (*t2>=0) && (*t2<=1))
        return(PARM_LINE_INTERSECT_IN_SEGMENT);
    else
        return(PARM_LINE_INTERSECT_OUT_SEGMENT);
}

int
Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2, POINT2D_PTR pt) {
    float t1, t2, det_p1p2 = (p1->v.x*p2->v.y - p1->v.y*p2->v.x);
    if (fabs(det_p1p2) <= EPSILON_E5) return(PARM_LINE_NO_INTERSECT);
    t1 = (p2->v.x*(p1->p0.y - p2->p0.y) - p2->v.y*(p1->p0.x - p2->p0.x))/det_p1p2;
    t2 = (p1->v.x*(p1->p0.y - p2->p0.y) - p1->v.y*(p1->p0.x - p2->p0.x)) /det_p1p2;
    pt->x = p1->p0.x + p1->v.x*t1;
    pt->y = p1->p0.y + p1->v.y*t1;
    if ((t1>=0) && (t1<=1) && (t2>=0) && (t2<=1))
        return(PARM_LINE_INTERSECT_IN_SEGMENT);
    else
        return(PARM_LINE_INTERSECT_OUT_SEGMENT);
}
/******************************************************************************/


/** 根据指定的点以及它们之间的向量,初始化一个3D参数化直线结构,在函数内部生成
 *  this initializes a parametric 3d line
 *  note that the function computes v=p_pend - p_init, thus when t=0 the line p=p0+v*t = p0
 *  and whne t=1, p=p1, this way the segement is traced out from p0 to p1 via 0<= t <= 1
 */
void
Init_Parm_Line3D(POINT3D_PTR p_init,POINT3D_PTR p_term, PARMLINE3D_PTR p) {
    // start point
    VECTOR3D_INIT(&(p->p0), p_init);
    // end point
    VECTOR3D_INIT(&(p->p1),p_term);
    // now compute direction vector from p0->p1
    VECTOR3D_Build(p_init, p_term, &(p->v));
}
/** 计算一条参数化直线在参数t处的值,并将返回结果存储在指定的中.
 *  this function computers the value of the sent parametric line at the value of t
 *
 *  @param p  <#p description#>
 *  @param t  <#t description#>
 *  @param pt <#pt description#>
 */
void
Compute_Parm_Line3D(PARMLINE3D_PTR p, float t, POINT3D_PTR pt) {
    pt->x = p->p0.x + p->v.x*t;
    pt->y = p->p0.y + p->v.y*t;
    pt->z = p->p0.z + p->v.z*t;
}

/******************************************************************************/

//3D平面函数
/** 使用指定的点和法线来初始化一个3D平面,另外,该函数可以对指定的法线进行归一化,使其长度为1.0
 *  要进行归一化,需要将normalize参数设置为True;否则,将它设置为FALSE.在很多光照和背面消除计算中.
 *  知道多边形或平面的法线长度为1.0会很有帮助.
 *  this function initializes a 3d plane
 *
 *  @param plane     <#plane description#>
 *  @param p0        <#p0 description#>
 *  @param normal    <#normal description#>
 *  @param normalize <#normalize description#>
 */
void
PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,VECTOR3D_PTR normal, int normalize) {
    // copy the point
    POINT3D_COPY(&plane->p0, p0);
    // if normalize is 1 then the normal is made into a unit vector
    if (!normalize)  VECTOR3D_COPY(&plane->n, normal);
    else  VECTOR3D_Normalize(normal,&plane->n);  // make normal into unit vector
}
/**  它判断指定点位于哪一个半空间,很多时候,需要判断某样东西位于平面的哪一边,该函数提供了这种功能。
 *   若指定点位于平面上返回0;   若位于正半空间中,返回一个正数;   若位于负半空间中,返回一个正数;
 *
 *  @param pt    <#pt description#>
 *  @param plane <#plane description#>
 *
 *  @return <#return value description#>
 */
float
Compute_Point_In_Plane3D(POINT3D_PTR pt, PLANE3D_PTR plane) {
    //检测点是在平面上 还是 正半空间 还是 负半空间
    //test if the point in on the plane, in the positive halfspace or negative halfspace
    float hs = plane->n.x*(pt->x - plane->p0.x) +
               plane->n.y*(pt->y - plane->p0.y) +
               plane->n.z*(pt->z - plane->p0.z);
    return hs;
}


/** 函数计算一条3D参数化直线与一个3D平面的交点,将交点处的参数值存储到 t中,并将交点存储到pt中。
 *  然而在使用这些数据之前,需要测试函数的返回值,以确定是否存在交点。
 *
 *  这个函数计算参数化直线与平面的交点
 *  计算交点时,该函数将参数化直线视为无穷长
 *  如果交点在线段pline上,则参数t的值将位于范围[0,1]内,
 *  另外,如果不相交,该函数返回0,如果交点在线段上,返回1
 *  如果返回交点不在线段上,返回2 ; 如果线段位于平面上,返回 3
 */
int
Intersect_Parm_Line3D_Plane3D(PARMLINE3D_PTR pline,PLANE3D_PTR plane,
                              float *t, POINT3D_PTR pt) {
    // 首先判断线段和平面是否平行
    // 若是,则它们不可能相交,除非线段位于平面上!
    float plane_dot_line = VECTOR3D_Dot(&pline->v, &plane->n);
    
    if (fabs(plane_dot_line) <= EPSILON_E5) {
        //线段与平面平行 ,它是否在平面上？
        if (fabs(Compute_Point_In_Plane3D(&pline->p0, plane)) <= EPSILON_E5)
            return(PARM_LINE_INTERSECT_EVERYWHERE);
        else
            return(PARM_LINE_NO_INTERSECT);
    }
    
    *t = -(plane->n.x*pline->p0.x +
           plane->n.y*pline->p0.y +
           plane->n.z*pline->p0.z -
           plane->n.x*plane->p0.x -
           plane->n.y*plane->p0.y -
           plane->n.z*plane->p0.z) / (plane_dot_line);
    // now plug t into the parametric line and solve for x,y,z
    pt->x = pline->p0.x + pline->v.x*(*t);
    pt->y = pline->p0.y + pline->v.y*(*t);
    pt->z = pline->p0.z + pline->v.z*(*t);
    // test interval of t [0,1]
    if (*t>=0.0 && *t<=1.0)
        return(PARM_LINE_INTERSECT_IN_SEGMENT );
    else
        return(PARM_LINE_INTERSECT_OUT_SEGMENT);
}
/******************************************************************************/

//定点数函数
/**
 *  this function computes the product fp_prod = fp1*fp2
 *  using 64 bit math, so as not to loose precission
 *  @param fp1 <#fp1 description#>
 *  @param fp2 <#fp2 description#>
 *
 *  @return <#return value description#>
 */
/*
FIXP16
FIXP16_MUL(FIXP16 fp1,FIXP16 fp2) {
    FIXP16 fp_prod; // return the product
    
    _asm {
        mov eax, fp1      // move into eax fp2
        imul fp2          // multiply fp1*fp2
        shrd eax, edx, 16 // result is in 32:32 format
                          // residing at edx:eax
                          // shift it into eax alone 16:16
                          // result is sitting in eax
    }
   
}
 */
/**
 *  this function computes the quotient fp1/fp2 using
 *  64 bit math, so as not to loose precision
 *  @param fp1 <#fp1 description#>
 *  @param fp2 <#fp2 description#>
 *
 *  @return <#return value description#>
 */
/*
FIXP16
FIXP16_DIV(FIXP16 fp1,FIXP16 fp2) {
    _asm {
        mov eax, fp1      // move dividend into eax
        cdq               // sign extend it to edx:eax
        shld edx, eax, 16 // now shift 16:16 into position in edx
        sal eax, 16       // and shift eax into position since the
                          // shld didn't move it -- DUMB! uPC
        idiv fp2          // do the divide
                          // result is sitting in eax
    }
}
 */

/******************************************************************************/


























