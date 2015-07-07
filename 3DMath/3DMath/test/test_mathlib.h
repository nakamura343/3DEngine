//
//  test_mathlib.h
//  3DMath
//
//  Created by TTc on 15/7/6.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#ifndef ___DMath__test_mathlib__
#define ___DMath__test_mathlib__

#include <stdio.h>
#include "T3DLib1.h"
#include "T3DLib4.h"



void test_MATRIX3X3();

void test_Fast_sin();

void test_Fast_cos();


void test_Init_Parm_Line2d();

void test_Compute_Parm_Line2D();

void test_Intersect_Parm_Lines2D();

void test_Init_Parm_Line3D();

void test_Compute_Parm_Line3D();

void test_PLANE3D_Init();

void test_Compute_Point_In_Plane3D();

void test_Intersect_Parm_Line3D_Plane3D();

void test_VECTOR3D_Theta_To_QUAT();

void test_EulerZYX_To_QUAT();

void test_QUAT_To_VECTOR3D_Theta();

void test_QUAT_Add();

void test_QUAT_Sub();

void test_QUAT_Conjugate();

void test_QUAT_Scale();

void test_QUAT_Normalize();

void test_QUAT_Unit_Inverse();

void test_QUAT_Inverse();

void test_QUAT_Mul();

void test_QUAT_Triple_Product();

#endif /* defined(___DMath__test_mathlib__) */
