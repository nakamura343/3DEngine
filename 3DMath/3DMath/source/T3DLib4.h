//
//  T3DLib.h
//  3DMathLib
//
//  Created by TTc on 15/7/2.
//  Copyright (c) 2015年 TTc. All rights reserved.
//
//
#ifndef ___DMathLib__T3DLib__
#define ___DMathLib__T3DLib__

//在计算机数学中 点与向量  几乎等效(至少是数据上等效)  LHS(左手坐标系)

#include <stdio.h>
#include <string.h>

/******************************************************************************/
// DEFINES & CONSTANTS //
// defines for small numbers
/*
 数学常量    16.16格式
 */


//与定点数运算相关的常量
#define FIXP16_SHIFT        16
#define FIXP16_MAG          65536
#define FIXP16_DP_MASK      0x0000ffff
#define FIXP16_WP_MASK      0xffff0000
#define FIXP16_ROUND_UP     0x00008000


//针对非常小的数的常量  Epsilon 常量有助于很小的浮点数进行数学比较
// eg: 执行浮点数运算时,经常需要检测一个浮点数是否等于0.0,然而经过几次浮点运算后,精度将降低,很少会出现0.0
//情况,因此我们只需要 检测是否接近于0.0
//if(fabs(x) < 0.00001f) { //足够小于0 }
#define EPSILON_E3          (float)(1E-3)
#define EPSILON_E4          (float)(1E-4)
#define EPSILON_E5          (float)(1E-5)
#define EPSILON_E6          (float)(1E-6)

//用于参数化直线交点的常量
#define PARM_LINE_NO_INTERSECT              0 //不相交
#define PARM_LINE_INTERSECT_IN_SEGMENT      1 //交点在线段上
#define PARM_LINE_INTERSECT_OUT_SEGMENT     2 //交点不在线段上
#define PARM_LINE_INTERSECT_EVERYWHERE      3 //代表 点位于平面内

/******************************************************************************/
/******************************************************************************/

//不包含w分量的 2d向量 和 2d点
typedef struct VECTOR2D_TYP {
    union{
        float M[2];     //数组存储方式
        struct{         //独立变量存储方式
            float x,y;
        };
    };
}VECTOR2D, POINT2D, *VECTOR2D_PTR , *POINT2D_PTR;

//不包含w分量的 3d向量 和 3d点
typedef struct VECTOR3D_TYP {
    union{
        float M[3];
        struct {
            float x,y,z;
        };
    };
}VECTOR3D, *VECTOR3D_PTR, POINT3D , *POINT3D_PTR;

//包含w分量的 4d向量 和 4d点
typedef struct VECTOR4D_TYP {
    union{
        float M[4];
        struct{
            float x,y,z,w;
        };
    };
}VECTOR4D, *VECTOR4D_PTR, POINT4D , *POINT4D_PTR;


/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/
//参数化直线

//2d 参数化直线  2D parametric line
typedef struct PARMLINE2D_TYP {  // |v| = |p0->p1 |
    POINT2D     p0; //参数化直线的起点
    POINT2D     p1; //参数化直线的终点
    VECTOR2D    v;  //线段的方向向量
}PARMLINE2D, *PARMLINE2D_PTR;


//3D parametric line
typedef struct PARMLINE3D_TYP {
    POINT3D     p0;
    POINT3D     p1;
    VECTOR3D    v;
    
}PARMLINE3D, *PARMLINE3D_PTR;


/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/
/*                  3D平面
 虽然在大多数时候,需要考虑3D平面,但通常是在多边形语境下考虑.由于以下几种原因,创建一个3d平面
 类型是个不错的主意:
 以便能够执行空间算法,简化裁剪算法的编写,有助于编写3d直线和3d平面之间的碰撞检测算法。
 表示3d平面的方法有很多,但归根结底它们表示的都是相同的东西.
 我将使用"点－法线" 形式(而不是显示形式) 来表达3d平面.然而这并不意味着不能以下形势存储平面
 
 a * x + b * y + c * z + d = 0
 
 这是一个函数,将它转换为 点－法线形式,方面存储,方面很多计算
 */
// 3D plane
typedef struct  PLANE3D_TYP {
    POINT3D     p0; //平面上的点
    VECTOR3D    n;  //平面的法线(不必是单位向量)
}PLANE3D, *PLANE3D_PTR;
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/******************************************************************************/

/*
 矩阵
 */

//4 X 4
typedef struct MATRIX4X4_TYP {
    union {
        float M[4][4];      //以数组形式存储
        struct {            //按照先行后列的 顺序 以独立变量的方式存储
            float M00,M01,M02,M03;
            float M10,M11,M12,M13;
            float M20,M21,M22,M23;
            float M30,M31,M32,M33;
        };
    };
}MATRIX4X4,  *MATRIX4X4_PTR;

//4 X 3
typedef struct MATRIX4X3_TYP {
    union {
        float M[4][3];
        struct{
            float M00,M01,M02;
            float M10,M11,M12;
            float M20,M21,M22;
            float M30,M31,M32;
        };
    };
}MATRIX4X3, *MATRIX4X3_PTR;

//1 X 4
typedef struct MATRIX1X4_TYP {
    union {
        float M[4];
        struct{
            float M00,M01,M02,M03;
        };
    };
}MATRIX1X4, *MATRIX1X4_PTR;

//3 X 3
typedef struct  MATRIX3X3_TYP {
    union{
        float M[3][3];
        struct{
            float M00,M01,M02;
            float M10,M11,M12;
            float M20,M21,M22;
        };
    };
}MATRIX3X3, *MATRIX3X3_PTR;

//1 X 3
typedef struct MATRIX1X3_TYP {
    union{
        float M[3];
        struct{
            float M00,M01,M02;
        };
    };
}MATRIX1X3, *MATRIX1X3_PTR;

// 3 X 2
typedef struct MATRIX3X2_TYP {
    union{
        float M[3][2];
        struct{
            float M00,M01;
            float M10,M11;
            float M20,M21;
        };
    };
}MATRIX3X2, *MATRIX3X2_PTR;

//2 X 2
typedef struct MATRIX2X2_TYP{
    union{
        float  M[2][2];
        struct{
            float M00,M01;
            float M10,M11;
        };
    };
}MATRIX2X2, *MATRIX2X2_PTR;

//1 X 2
typedef struct MATRIX1X2_TYP{
    union{
        float M[2];
        struct{
            float M00,M01;
        };
    };
}MATRIX1X2, *MATRIX1X2_PTR;


/******************************************************************************/
/******************************************************************************/
/*     四元数
 四元数由4个分量组成,实部为q0,向量部分为虚部.另外我们将在旋转和相机语境下使用四元数,所有用实数
 q0(或w)和向量部分<x,y,z>来表示四元数有很多的优点。
 
 这样可以使用其他VECTOR3D函数来对四元数执行操作--可能由于他们的数据结构类似
 */
//四元数        4d quaternion
typedef struct QUAT_TYP {
    union{
        float M[4];  //按w,x,y,z 的顺序以 数组方式存储
        struct {     //向量形式
            float       q0;           //实部
            VECTOR3D    qv;           //虚部
        };
        struct {     //显式
            float w,x,y,z;
        };
    };
}QUAT, *QUAT_PTR;


/******************************************************************************/
/******************************************************************************/
/*  2D极坐标
 
 构建炮塔模型使用极坐标系将最合适,因为需要相对于坐标轴的角度方向,因此，数学引擎支持所有的2d和3d
 角坐标系以及这些坐标系之间的变换。 支持的坐标系 包括 2d极坐标系，3d柱面坐标系 ，3d球面坐标系
 所有角度都以弧度为单位，所有转换都是针对 RHS(右手坐标系)的,如果你使用左手坐标系,一定要小心
 
 如果不考虑坐标系的右手性,得到的结果可能是正确结果求负.
 */

//2D极坐标  // 2D polar coordinates
typedef struct POLAR2D_TYP {
    float   r;        //半径
    float   theta;    //角度
}POLAR2D, *POLAR2D_PTR;



//3d 柱面坐标    3D cylindrical coordinates
typedef struct CYLINDRICAL3D_TYP {
    float   r;        //半径
    float   theta;    //与z轴的夹角
    float   z;        //坐标
    
}CYLINDRICAL3D, *CYLINDRICAL3D_PTR;

//3d 球面坐标系食所有角坐标系中最复杂的。 3D spherical coordinates
typedef struct SPHERICAL3D_TYP {
    float   p;         //到原点的距离
    float   theta;     //线段0->p 和 正z轴 之间的夹角
    float   phi;       //线段0->p 和 在x－y平面上的投影 与正x轴 之间的夹角
}SPHERICAL3D, *SPHERICAL3D_PTR;

/******************************************************************************/
/******************************************************************************/

/*   定点数 ->非常难
 定点数用于以有限的精度表示浮点数的整数,以避免进行浮点处理. 现在对它们的需求已经不再像以前
 那么强烈,因为cpu已经可以像定点数一样 迅速地处理执行 浮点数运算。
 然而，有时候仍需要使用定点数, eg: 无法使用浮点的内部循环，在你只想使用整数时或者在小型手持设备
 上时, 以下使用16.16格式, 即其中使用16bit 表示整数部分,余下的16bit表示小数部分.
 */
//定点数类型
typedef int FIXP16;

typedef int *FIXP16_PTR;


/******************************************************************************/
/******************************************************************************/





//单位矩阵
//4 X 4 单位矩阵
const MATRIX4X4 IMAT_4X4 = {1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1};

//4 x 3 单位矩阵  (从数学上说是不正确的,本书后面将使用这种矩阵,并假设最后一列为[0 0 0 1] t) 只有 nXn 方阵才有单位矩阵
//  我们使用4x3矩阵，并假设最后一列为[0 0 0 1] t ,因此，可以定义一个 4 X 3单位矩阵，但在执行
//数学运算时 需要记住它有一个额外的列.
const MATRIX4X3 IMAT_4X3 = {1,0,0,
    0,1,0,
    0,0,1,
    0,0,0};

//3 X 3 单位矩阵
const MATRIX3X3 IMAT_3X3 = {1,0,0,
    0,1,0,
    0,0,1};

//2 X 2 单位矩阵
const MATRIX2X2 IMAT_2X2 = {1,0,
    0,1};
/******************************************************************************/
/******************************************************************************/



//点和向量宏

//向量归零宏
inline void
VECTOR2D_ZERO(VECTOR2D_PTR v) {
    (v)->x = (v)->y = 0.0;
}

inline void
VECTOR3D_ZERO(VECTOR3D_PTR v) {
    (v)->x = (v)->y = (v)->z = 0.0;
}

//向量宏,4D向量的w被设置为1
inline void
VECTOR4D_ZERO(VECTOR4D_PTR v) {
    (v)->x = (v)->y = (v)->z = 0.0;
    (v)->w = 1.0;
}

//使用分量初始化向量的宏
inline void
VECTOR2D_INITXY(VECTOR2D_PTR v,float x, float y) {
    (v)->x = x;
    (v)->y = y;
}

//使用分量初始化向量的宏
inline void
VECTOR3D_INITXYZ(VECTOR3D_PTR v,float x, float y,float z) {
    (v)->x = x;
    (v)->y = y;
    (v)->z = z;
}

//使用分量初始化向量的宏
inline void
VECTOR4D_INITXYZ(VECTOR4D_PTR v,float x, float y,float z) {
    (v)->x = x;
    (v)->y = y;
    (v)->z = z;
    (v)->w = 1.0;
}

//使用另一个向量来初始化向量的宏
inline void
VECTOR2D_INIT(VECTOR2D_PTR vdst,VECTOR2D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
}

inline void
VECTOR3D_INIT(VECTOR3D_PTR vdst,VECTOR3D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void
VECTOR4D_INIT(VECTOR4D_PTR vdst,VECTOR4D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
    (vdst)->w = (vsrc)->w;
}

//复制向量的宏
inline void
VECTOR2D_COPY(VECTOR2D_PTR vdst,VECTOR2D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
}


inline void
VECTOR3D_COPY(VECTOR3D_PTR vdst,VECTOR3D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void
VECTOR4D_COPY(VECTOR4D_PTR vdst,VECTOR4D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
    (vdst)->w = (vsrc)->w;
}

//初始化点的宏
inline void
POINT2D_INIT(POINT2D_PTR vdst,POINT2D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
}

inline void
POINT3D_INIT(POINT3D_PTR vdst,POINT3D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void
POINT4D_INIT(POINT4D_PTR vdst,POINT4D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
    (vdst)->w = (vsrc)->w;
}

//复制点的宏
inline void
POINT2D_COPY(POINT2D_PTR vdst,POINT2D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
}

inline void
POINT3D_COPY(POINT3D_PTR vdst,POINT3D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void
POINT4D_COPY(POINT4D_PTR vdst,POINT4D_PTR vsrc) {
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
    (vdst)->w = (vsrc)->w;
}

//矩阵宏

//清空矩阵的宏
#define MAT_ZERO_2X2(m)  {memset((void*)(m),0,sizeof(MATRIX2X2));}
#define MAT_ZERO_3X3(m)  {memset((void*)(m),0,sizeof(MATRIX3X3));}
#define MAT_ZERO_4X4(m)  {memset((void*)(m),0,sizeof(MATRIX4X4));}
#define MAT_ZERO_4X3(m)  {memset((void*)(m),0,sizeof(MATRIX4X3));}


//设置单位矩阵的宏
#define MAT_IDENTITY_2X2(m)  {memcpy((void*)(m),(void*)&IMAT_2X2,sizeof(MATRIX2X2));}
#define MAT_IDENTITY_3X3(m)  {memcpy((void*)(m),(void*)&IMAT_3X3,sizeof(MATRIX3X3));}
#define MAT_IDENTITY_4X4(m)  {memcpy((void*)(m),(void*)&IMAT_4X4,sizeof(MATRIX4X4));}
#define MAT_IDENTITY_4X3(m)  {memcpy((void*)(m),(void*)&IMAT_4X3,sizeof(MATRIX4X3));}


//复制矩阵的宏
#define MAT_COPY_2X2(src_mat,dest_mat) {memcpy((void*)(dest_mat),\
(void*)(src_mat),sizeof(MATRIX2X2)); }

#define MAT_COPY_3X3(src_mat,dest_mat) {memcpy((void*)(dest_mat),\
(void*)(src_mat),sizeof(MATRIX3X3)); }

#define MAT_COPY_4X4(src_mat,dest_mat) {memcpy((void*)(dest_mat),\
(void*)(src_mat),sizeof(MATRIX4X4)); }

#define MAT_COPY_4X3(src_mat,dest_mat) {memcpy((void*)(dest_mat),\
(void*)(src_mat),sizeof(MATRIX4X3)); }



//矩阵进行转置的宏   matrix transposing macros
inline void
MAT_TRANSPOSE_3X3(MATRIX3X3_PTR m) {
    MATRIX3X3 mt;
    mt.M00 = m->M00;mt.M01 = m->M10;mt.M02 = m->M20;
    mt.M10 = m->M01;mt.M11 = m->M11;mt.M12 = m->M21;
    mt.M20 = m->M02;mt.M21 = m->M12;mt.M22 = m->M22;
    memcpy((void *)m, (void *)&mt, sizeof(MATRIX3X3));
}


inline void
MAT_TRANSPOSE_4X4(MATRIX4X4_PTR m) {
    MATRIX4X4 mt;
    mt.M00 = m->M00;mt.M01 = m->M10;mt.M02 = m->M20; mt.M03 = m->M30;
    mt.M10 = m->M01;mt.M11 = m->M11;mt.M12 = m->M21; mt.M13 = m->M31;
    mt.M20 = m->M02;mt.M21 = m->M12;mt.M22 = m->M22; mt.M23 = m->M32;
    mt.M30 = m->M03;mt.M31 = m->M13;mt.M32 = m->M23; mt.M33 = m->M33;
    memcpy((void *)m, (void *)&mt, sizeof(MATRIX3X3));
}

inline void
MAT_TRANSPOSE_3X3_PTR(MATRIX3X3_PTR m,MATRIX3X3_PTR mt) {
    mt->M00 = m->M00;mt->M01 = m->M10;mt->M02 = m->M20;
    mt->M10 = m->M01;mt->M11 = m->M11;mt->M12 = m->M21;
    mt->M20 = m->M02;mt->M21 = m->M12;mt->M22 = m->M22;
}

inline void
MAT_TRANSPOSE_4X4_PTR(MATRIX4X4_PTR m,MATRIX4X4_PTR mt) {
    mt->M00 = m->M00;mt->M01 = m->M10;mt->M02 = m->M20; mt->M03 = m->M30;
    mt->M10 = m->M01;mt->M11 = m->M11;mt->M12 = m->M21; mt->M13 = m->M31;
    mt->M20 = m->M02;mt->M21 = m->M12;mt->M22 = m->M22; mt->M23 = m->M32;
    mt->M30 = m->M03;mt->M31 = m->M13;mt->M32 = m->M23; mt->M33 = m->M33;
}

//小型内联函数,可以将其转化为宏

//矩阵和向量列 互换宏  (所有的矩阵函数都使用 指针)
inline void
MAT_COLUMN_SWAP_2X2(MATRIX2X2_PTR m, int c,MATRIX1X2_PTR v) {
    m->M[0][c] = v->M[0];
    m->M[1][c] = v->M[1];
}

inline void
MAT_COLUMN_SWAP_3X3(MATRIX3X3_PTR m, int c,MATRIX1X3_PTR v) {
    m->M[0][c] = v->M[0];
    m->M[1][c] = v->M[1];
    m->M[2][c] = v->M[2];
}

inline void
MAT_COLUMN_SWAP_4X4(MATRIX3X3_PTR m, int c,MATRIX1X4_PTR v) {
    m->M[0][c] = v->M[0];
    m->M[1][c] = v->M[1];
    m->M[2][c] = v->M[2];
    m->M[3][c] = v->M[3];
}

inline void
MAT_COLUMN_SWAP_4X3(MATRIX4X3_PTR m, int c,MATRIX1X4_PTR v) {
    m->M[0][c] = v->M[0];
    m->M[1][c] = v->M[1];
    m->M[2][c] = v->M[2];
    m->M[3][c] = v->M[3];
}

//四元数宏  根据需要的格式初始化四元数的宏  quaternion macros
inline void
QUAT_ZERO(QUAT_PTR q) {
    (q)->x = (q)->y = (q)->z = (q)->w = 0.0;
}

inline void
QUAT_INITWXYZ(QUAT_PTR q,float w,float x,float y,float z) {
    (q)->w = w;(q)->x = x; (q)->y = y; (q)->z = z;
}

inline void
QUAT_INIT_VECTOR3D(QUAT_PTR q,VECTOR3D_PTR v) {
    (q)->w = 0; (q)->x = v->x; (q)->y = v->y; (q)->z = v->z;
}

inline void
QUAT_INIT(QUAT_PTR qdst, QUAT_PTR qsrc) {
    (qdst)->w = (qsrc)->w;(qdst)->x = (qsrc)->x;
    (qdst)->y = (qsrc)->y;(qdst)->z = (qsrc)->z;
}

inline void
QUAT_COPY(QUAT_PTR qdst,QUAT_PTR qsrc) {
    (qdst)->w = (qsrc)->w;(qdst)->x = (qsrc)->x;
    (qdst)->y = (qsrc)->y;(qdst)->z = (qsrc)->z;
}

//定点数宏
#define FIXP16_WP(fp) ((fp) >> FIXP16_SHIFT)
#define FIXP16_DP(fp) ((fp) && FIXP16_DP_MASK)
//将整数和浮点数转换为16.16格式的定点数
#define INT_TO_FIXP16(i) ((i) << FIXP16_SHIFT)
#define FLOAT_TO_FIXP16(f) (((float)(f) * (float)FIXP16_MAG + 0.5))
//将定点数转换为浮点数
#define FIXP16_TO_FLOAT(fp) (((float)fp)/FIXP16_MAG)

/******************************************************************************/
/******************************************************************************/

/*
 函数原型   FUNC PROTOTYPES
 */

//通用三角函数
float Fast_Sin(float theta);

float Fast_Cos(float theta);

//极坐标,柱面坐标,球面坐标
void POLAR2D_To_POINT2D(POLAR2D_PTR polar, POINT2D_PTR rect);

void POLAR2D_To_RectXY(POLAR2D_PTR polar, float *x, float *y);

void POINT2D_To_POLAR2D(POINT2D_PTR rect,POLAR2D_PTR polar);

void POINT2D_To_PolarRTh(POINT2D_PTR rect,float *r , float *theta);

//柱面坐标 函数
void CYLINDRICAL3D_To_POINT3D(CYLINDRICAL3D_PTR cyl,POINT3D_PTR rect);

void CYLINDRICAL3D_To_RectXYZ(CYLINDRICAL3D_PTR cyl,float *x,float *y,float *z);

void POINT3D_To_CYLINDRICAL3D(POINT3D_PTR rect,CYLINDRICAL3D_PTR cyl);

void POINT3D_To_CylindricalRThZ(POINT3D_PTR rect,float *r ,float *theta, float *z);

//球面坐标 函数
void SPHERICAL3D_To_POINT3D(SPHERICAL3D_PTR sph,POINT3D_PTR rect);

void SPHERICAL3D_To_RectXYZ(SPHERICAL3D_PTR sph, float *x,float *y,float *z);

void POINT3D_To_SPHERICAL3D(POINT3D_PTR rect,SPHERICAL3D_PTR sph);

void POINT3D_To_SphericalPThPh(POINT3D_PTR rect,float *p ,float *theta, float *phi);


//2D向量函数
void VECTOR2D_ADD(VECTOR2D_PTR va, VECTOR2D_PTR vb, VECTOR2D_PTR vsum);

VECTOR2D VECTOR2D_ADD(VECTOR2D_PTR va , VECTOR2D_PTR vb);

void VECTOR2D_Sub(VECTOR2D_PTR va, VECTOR2D_PTR vb, VECTOR2D_PTR vdiff);

VECTOR2D VECTOR2D_Sub(VECTOR2D_PTR va, VECTOR2D_PTR vb) ;

void VECTOR2D_Scale(float k,VECTOR2D_PTR va);

void VECTOR2D_Scale(float k,VECTOR2D_PTR va,VECTOR2D_PTR vscaled);

float VECTOR2D_Dot(VECTOR2D_PTR va,VECTOR2D_PTR vb);

float VECTOR2D_Length(VECTOR2D_PTR va);

float VECTOR2D_Length_Fast(VECTOR2D_PTR va);

void VECTOR2D_Normalize(VECTOR2D_PTR va);

void VECTOR2D_Normalize(VECTOR2D_PTR va,VECTOR2D_PTR vn);

void VECTOR2D_Build(VECTOR2D_PTR init ,VECTOR2D_PTR term , VECTOR2D_PTR result);

float VECTOR2D_CosTh(VECTOR2D_PTR va, VECTOR2D_PTR vb);

void VECTOR2D_Print(VECTOR2D_PTR va, char *name);


//3D向量函数
void VECTOR3D_ADD(VECTOR3D_PTR va, VECTOR3D_PTR vb, VECTOR3D_PTR vsum);

VECTOR3D VECTOR3D_ADD(VECTOR3D_PTR va , VECTOR3D_PTR vb);

void VECTOR3D_Sub(VECTOR3D_PTR va, VECTOR3D_PTR vb, VECTOR3D_PTR vdiff);

VECTOR3D VECTOR3D_Sub(VECTOR3D_PTR va, VECTOR3D_PTR vb);

void VECTOR3D_Scale(float k,VECTOR3D_PTR va);

void VECTOR3D_Scale(float k,VECTOR3D_PTR va,VECTOR3D_PTR vscaled);

float VECTOR3D_Dot(VECTOR3D_PTR va,VECTOR3D_PTR vb);

void VECTOR3D_Cross(VECTOR3D_PTR va,VECTOR3D_PTR vb,VECTOR3D_PTR vn);

VECTOR3D VECTOR3D_Cross(VECTOR3D_PTR va,VECTOR3D_PTR vb);

float VECTOR3D_Length(VECTOR3D_PTR va);

float VECTOR3D_Length_Fast(VECTOR3D_PTR va);

void VECTOR3D_Normalize(VECTOR3D_PTR va);

void VECTOR3D_Normalize(VECTOR3D_PTR va,VECTOR3D_PTR vn);

void VECTOR3D_Build(VECTOR3D_PTR init ,VECTOR3D_PTR term , VECTOR3D_PTR result);

float VECTOR3D_CosTh(VECTOR3D_PTR va, VECTOR3D_PTR vb);

void VECTOR3D_Print(VECTOR3D_PTR va, char *name);


//4D向量函数
void VECTOR4D_ADD(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum);

VECTOR4D VECTOR4D_ADD(VECTOR4D_PTR va , VECTOR4D_PTR vb);

void VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vdiff);

VECTOR4D VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb);

void VECTOR4D_Scale(float k,VECTOR4D_PTR va);

void VECTOR4D_Scale(float k,VECTOR4D_PTR va,VECTOR4D_PTR vscaled);

float VECTOR4D_Dot(VECTOR4D_PTR va,VECTOR4D_PTR vb);

void VECTOR4D_Cross(VECTOR4D_PTR va,VECTOR4D_PTR vb,VECTOR4D_PTR vn);

VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va,VECTOR4D_PTR vb);

float VECTOR4D_Length(VECTOR4D_PTR va);

float VECTOR4D_Length_Fast(VECTOR4D_PTR va);

void VECTOR4D_Normalize(VECTOR4D_PTR va);

void VECTOR4D_Normalize(VECTOR4D_PTR va,VECTOR4D_PTR vn);

void VECTOR4D_Build(VECTOR4D_PTR init ,VECTOR4D_PTR term , VECTOR4D_PTR result);

float VECTOR4D_CosTh(VECTOR4D_PTR va, VECTOR4D_PTR vb);

void VECTOR4D_Print(VECTOR4D_PTR va, char *name);

/******************************************************************************/
// 2x2 matrix functions (note there others in T3DLib1.cpp|h)
void Mat_Init_2X2(MATRIX2X2_PTR ma,float m00,float m01,float m10,float m11);

void Print_Mat_2X2(MATRIX2X2_PTR ma, char *name);

float Mat_Det_2X2(MATRIX2X2_PTR m);

void Mat_Add_2X2(MATRIX2X2_PTR ma,MATRIX2X2_PTR mb,MATRIX2X2_PTR msum);

void Mat_Mul_2X2(MATRIX2X2_PTR ma,MATRIX2X2_PTR mb,MATRIX2X2_PTR mprod);

int Mat_Inverse_2X2(MATRIX2X2_PTR ma,MATRIX2X2_PTR mi);
//
int Solve_2X2_System(MATRIX2X2_PTR A,MATRIX1X2_PTR X,MATRIX1X2_PTR B);


//3X3 matrix functions
int MAT_Mul_1X2_3X2(MATRIX1X2_PTR ma, MATRIX3X2_PTR mb,MATRIX1X2_PTR mprod);

int Mat_Mul_1X3_3X3(MATRIX1X3_PTR ma, MATRIX3X3_PTR mb,MATRIX1X3_PTR mprod);

int Mat_Mul_3X3(MATRIX3X3_PTR ma,MATRIX3X3_PTR mb,MATRIX3X3_PTR mprod);


inline int Mat_Init_3X2(MATRIX3X2_PTR ma,float m00, float m01,float m10, float m11,
                        float m20, float m21);


void Mat_Add_3X3(MATRIX3X3_PTR ma, MATRIX3X3_PTR mb, MATRIX3X3_PTR msum);

void Mat_Mul_VECTOR3D_3X3(VECTOR3D_PTR  va, MATRIX3X3_PTR mb,VECTOR3D_PTR  vprod);

void Mat_Init_3X3(MATRIX3X3_PTR ma,float m00, float m01, float m02,float m10, float m11, float m12,
                  float m20, float m21, float m22);

int Mat_Inverse_3X3(MATRIX3X3_PTR m, MATRIX3X3_PTR mi);

float Mat_Det_3X3(MATRIX3X3_PTR m);

int Solve_3X3_System(MATRIX3X3_PTR A, MATRIX1X3_PTR X, MATRIX1X3_PTR B);

//4x4 matrix functions
//ma + mb 将俩个矩阵相加,并将结果存储懂msum
void Mat_Add_4X4(MATRIX4X4_PTR ma,MATRIX4X4_PTR mb,MATRIX4X4_PTR msum);
//ma * mb 将俩矩阵相乘,并将结果存储在mprod
void Mat_Mul_4X4(MATRIX4X4_PTR ma,MATRIX4X4_PTR mb,MATRIX4X4_PTR mprod);
//将一个1X4的行向量(4D点)于一个4x4矩阵相乘
void Mat_Mul_1X4_4X4(MATRIX1X4_PTR ma,MATRIX4X4_PTR mb,MATRIX1X4_PTR mprod);
//将一个3d向量(1X3矩阵)与一个4x4矩阵相乘
void Mat_Mul_VECTOR3D_4X4(VECTOR3D_PTR va,MATRIX4X4_PTR mb, VECTOR3D_PTR mprod);
//将一个3d向量(1X3矩阵)与一个4x3矩阵相乘 (假设向量va中存在第四个元素=1.0)
void Mat_Mul_VECTOR3D_4X3(VECTOR3D_PTR  va,MATRIX4X3_PTR mb,VECTOR3D_PTR  vprod);
//将一个4d向量(1X4矩阵)与一个4x4矩阵相乘
void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR  va,MATRIX4X4_PTR mb,VECTOR4D_PTR  vprod);

void Mat_Mul_VECTOR4D_4X3(VECTOR4D_PTR  va,MATRIX4X4_PTR mb,VECTOR4D_PTR  vprod);
//使用传入的浮点值以先行后列的次序初始化矩阵ma.
void Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33);
//计算矩阵m的逆矩阵(如果存在的情况),若存在逆矩阵则返回1,否则返回0
int Mat_Inverse_4X4(MATRIX4X4_PTR m, MATRIX4X4_PTR mi);


//四元数函数 quaternions functions
/*
   四元素没有太多的物理意义,无法以图形方式进行检查。
   使用四元数会使很多操作简单化,eg:对于绕任意直线旋转的问题,如果使用四元数来求解将非常方便,
   而使用欧拉方程则非常复杂,四元数就是一个超复数w 为实部,另外3个为虚部
 */

//将两个四元数q1和q2 相加,结果存储到qsum中
void QUAT_Add(QUAT_PTR q1,QUAT_PTR q2,QUAT_PTR qsum);
////将两个四元数q1和q2 相减,结果存储到qdiff中
void QUAT_Sub(QUAT_PTR q1,QUAT_PTR q2,QUAT_PTR qdiff);
//计算四元数q的共轭,并将结果存储到qconj中
void QUAT_Conjugate(QUAT_PTR q,QUAT_PTR qconj);
//根据缩放因子scale对四元数q进行缩放,并将结果储存到qs中
void QUAT_Scale(QUAT_PTR q,float scale,QUAT_PTR qs);
//直接缩放四元素q,即修改q
void QUAT_Scale(QUAT_PTR q,float scale);
//返回四元数的范数,即长度
float QUAT_Norm(QUAT_PTR q);
//返回四元数q的范数的平方,即长度的平方
float QUAT_Norm2(QUAT_PTR q);
//函数将四元数q归一化,并将结果存储在qn中
void  QUAT_Normalize(QUAT_PTR q,QUAT_PTR qn);

void QUAT_Normalize(QUAT_PTR q);
//计算四元数q的逆,并将结果存储到qi中(q必须为单位四元数);因为该函数基于单位四元数的逆就是其共轭这一原理来求逆
void QUAT_Unit_Inverse(QUAT_PTR q, QUAT_PTR qi);

void QUAT_Unit_Inverse(QUAT_PTR q);
//计算非单位四元数q的逆
void QUAT_Inverse(QUAT_PTR q, QUAT_PTR qi);

void QUAT_Inverse(QUAT_PTR q);
//两个四元数相乘
void QUAT_Mul(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qprod);
//3个四元数相乘,并将结果存储到qprod中,该函数对  点旋转有很大用处
void QUAT_Triple_Product(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR q3,
                         QUAT_PTR qprod);
//根据方向向量v和角度theta创建一个旋转四元数。
void VECTOR3D_Theta_To_QUAT(QUAT_PTR q, VECTOR3D_PTR v, float theta);

void VECTOR4D_Theta_To_QUAT(QUAT_PTR q, VECTOR4D_PTR v, float theta);
//根据绕z,y,x 旋转的欧拉角创建一个旋转四元数.它是一种基本的相机转换.旋转顺序有6种(3的阶乘),zyz次序是最常见的
void EulerZYX_To_QUAT(QUAT_PTR q, float theta_z, float theta_y, float theta_x);
//将一个单位旋转四元数转换为一个单位3D向量和与一个绕该向量旋转的角度theta. VECTOR3D_Theta_To_QUAT()的反函数
void QUAT_To_VECTOR3D_Theta(QUAT_PTR q, VECTOR3D_PTR v, float *theta);

/******************************************************************************/

//2D参数化直线函数
void Init_Parm_Line2D(POINT2D_PTR p_init,POINT2D_PTR p_term,PARMLINE2D_PTR p);

void Compute_Parm_Line2D(PARMLINE2D_PTR p, float t, POINT2D_PTR pt);

int Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2,
                           float *t1, float *t2);

int Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2, POINT2D_PTR pt);

/******************************************************************************/

//3D参数化直线函数
void Init_Parm_Line3D(POINT3D_PTR p_init,POINT3D_PTR p_term, PARMLINE3D_PTR p);

void Compute_Parm_Line3D(PARMLINE3D_PTR p, float t, POINT3D_PTR pt);

/******************************************************************************/

//3D平面函数
void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,VECTOR3D_PTR normal, int normalize);

float Compute_Point_In_Plane3D(POINT3D_PTR pt, PLANE3D_PTR plane);

int Intersect_Parm_Line3D_Plane3D(PARMLINE3D_PTR pline,PLANE3D_PTR plane,
                                  float *t, POINT3D_PTR pt);

/******************************************************************************/
//定点数函数
// 2个定点数相乘
FIXP16 FIXP16_MUL(FIXP16 fp1,FIXP16 fp2);
// 2个定点数相除
FIXP16 FIXP16_DIV(FIXP16 fp1,FIXP16 fp2);


/******************************************************************************/


#endif /* defined(___DMathLib__T3DLib__) */
