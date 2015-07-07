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


void
test_Init_Parm_Line2d() {
    POINT2D p1 = {1,2}, p2 = {10,20};
    PARMLINE2D p;
    //创建一条从p1到p2的参数化直线
    Init_Parm_Line2D(&p1, &p2, &p);
}

void
test_Compute_Parm_Line2D(){
    POINT2D p1 = {1,2}, p2 = {10,20} ,pt;
    PARMLINE2D p;
    //创建一条从p1到p2的参数化直线
    Init_Parm_Line2D(&p1, &p2, &p);

    //计算t＝0.5对应的参数化直线上的点,并存储在pt中
    Compute_Parm_Line2D(&p, 0.5, &pt);
}


void
test_Intersect_Parm_Lines2D(){
    POINT2D p1 = {1,1} , p2 = {9,8};
    POINT2D p3 = {2,8} , p4 = {7,1};
    
    PARMLINE2D pl1,pl2;
    //创建一条从p1到p2的参数化直线
    Init_Parm_Line2D(&p1, &p2, &pl1);
    //创建一条从p3到p4的参数化直线
    Init_Parm_Line2D(&p3, &p4, &pl2);
    
    
    float t1 = 0.0, t2 = 0.0;
    //计算交点
    int intersection_type = Intersect_Parm_Lines2D(&pl1, &pl2, &t1,&t2);
    printf("intersection_type==%d",intersection_type);
}

void
test_Init_Parm_Line3D() {
    POINT3D p1 = {1,2,3}, p2 = {10,20,30};
    PARMLINE3D p;
    Init_Parm_Line3D(&p1, &p2, &p);
}

void
test_Compute_Parm_Line3D() {
    POINT3D p1 = {1,2,3}, p2 = {10,20,30},pt;
    PARMLINE3D p;
    
    Init_Parm_Line3D(&p1, &p2, &p);
    
    Compute_Parm_Line3D(&p, 0.5, &pt);
}


void
test_PLANE3D_Init() {
    VECTOR3D n = {1,1,1};
    POINT3D  p = {0,0,0};
    PLANE3D  plane;
    
    PLANE3D_Init(&plane, &p, &n, true);
}


void
test_Compute_Point_In_Plane3D() {
    VECTOR3D n = {1,1,1};
    POINT3D  p = {0,0,0};
    PLANE3D  plane;
    
    PLANE3D_Init(&plane, &p, &n, true);
    
    POINT3D p_test = {50,50,50};
    //检测点的位置
    float hs = Compute_Point_In_Plane3D(&p_test, &plane);
    printf("hs==%.2f",hs);
}

void
test_Intersect_Parm_Line3D_Plane3D() {
    POINT3D p1 = {5,5,-5}, p2 = {5,5,5},pt;
    PARMLINE3D pl;
    float t;
    //这条直线与z轴平行
    Init_Parm_Line3D(&p1, &p2, &pl);
    
    VECTOR3D n = {0,0,1};
    POINT3D p = {0,0,0};
    PLANE3D plane;
    //该平面与 x-y平面平行
    PLANE3D_Init(&plane, &p, &n, true);
    
    //计算交点;应为(5,5,0)
    int intersection_type = Intersect_Parm_Line3D_Plane3D(&pl, &plane, &t, &pt);
    printf("intersection_type==%d",intersection_type);

    
}

