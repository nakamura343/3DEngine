//
//  T3DLib5.cpp
//  3DMath
//
//  Created by TTc on 15/4/10.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#include "T3DLib5.h"
#include "T3DLib1.h"

#include <math.h>
//调用重置函数,以重置渲染列表,供动画的下一帧使用
void
Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list) {
    //这个函数初始化和重置传递进来的渲染列表
    //为将多边形插入到其中做好准备
    //这个版本的渲染列表由一个FACE4DV1 指针数组组成
    //每一帧都需要调用 这个函数
    //这个使用 num_polys 来跟踪渲染列表中包含的多边形数目
    //因此将其设置为0
    //如果需要使渲染列表更通用,需要采用更健壮的方案
    //并将其与多边形指针列表的关联切断
    rend_list->num_polys = 0;   //硬编码
}


/*  这个函数使用传递进来的矩阵;对渲染列表中局部顶点数组或变换后的顶点数组中所有的多边形顶点进行变换
 *
 */
void
Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,MATRIX4X4_PTR mt,int coord_select) {
    
    //应对哪个数组中的坐标进行变换
    switch(coord_select){
        case TRANSFORM_LOCAL_ONLY: {
            for (int poly = 0; poly < rend_list->num_polys; poly++) {
                //获得当前多边形
                POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
                //判断这个多边形是否有效？;当且仅当多边形处于活动状态并可见时才对其进行变换
                //在线框引擎中,“背面”的概念无关紧要
                if ((curr_poly==NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                    (curr_poly->state & POLY4DV1_STATE_CLIPPED ) ||
                    (curr_poly->state & POLY4DV1_STATE_BACKFACE) )
                    continue;
                //满足条件,进行变换
                for (int vertex = 0; vertex < 3; vertex++) {
                    //使用矩阵mt对顶点进行变换;用于暂时储存变换结果
                    POINT4D presult;
                    //对点进行变换
                    Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &presult);
                    //将结构存回去
                    VECTOR4D_COPY(&curr_poly->vlist[vertex], &presult);
                }
            }
        }
            break;
            
        case TRANSFORM_TRANS_ONLY:  {
            //对渲染列表中每个变换后的顶点进行变换;数组tvlist[]用于储存累积变换结果
            for (int poly = 0; poly < rend_list->num_polys; poly++) {
                //获得当前多边形
                POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
                //盘对
                if ((curr_poly==NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                    (curr_poly->state & POLY4DV1_STATE_CLIPPED ) ||
                    (curr_poly->state & POLY4DV1_STATE_BACKFACE) )
                    continue;
                
                for (int vertex = 0; vertex < 3; vertex++) {
                    POINT4D presult;
                    Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex], mt, &presult);
                    VECTOR4D_COPY(&curr_poly->tvlist[vertex], &presult);
                }
            }
        }
            break;
            
        case TRANSFORM_LOCAL_TO_TRANS: {
            //对渲染列表中的局部/模型顶点列表进行变换;并将结果存储到变换后的顶点列表中
            for (int poly = 0; poly < rend_list->num_polys; poly++) {
                //获得下一个多边形
                POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
                
                if ((curr_poly==NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                    (curr_poly->state & POLY4DV1_STATE_CLIPPED ) ||
                    (curr_poly->state & POLY4DV1_STATE_BACKFACE) )
                    continue;
                
                for (int vertex = 0; vertex < 3; vertex++) {
                    //使用矩阵mt对顶点进行变换
                    Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &curr_poly->tvlist[vertex]);
                }
            }
        }
            break;
        default: break;
    }
}
//
/* 这个函数使用传递进来的矩阵;对局部数组或变换后的数组中的所有顶点进行变换
 *
 *
 *
 */
void
Transform_OBJECT4DV1(OBJECT4DV1_PTR obj,MATRIX4X4_PTR mt, int coord_select, int transform_basis) {
    //判断对哪种坐标进行变换
    switch(coord_select) {
        case TRANSFORM_LOCAL_ONLY: {
            //对物体的每个局部/模型顶点坐标进行变换
            for (int vertex=0; vertex < obj->num_vertices; vertex++) {
                POINT4D presult;//用于暂时储存变换结果
                //对顶点进行变换
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex], mt, &presult);
                //将结果存回去
                VECTOR4D_COPY(&obj->vlist_local[vertex], &presult);
            }
        }
            break;
            
        case TRANSFORM_TRANS_ONLY: {
            //对物体的每个变换后的顶点进行变换,数组vlist_trans[]用于储存累积变换结果
            for (int vertex=0; vertex < obj->num_vertices; vertex++) {
                POINT4D presult;
                //对顶点进行变换
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex], mt, &presult);
                VECTOR4D_COPY(&obj->vlist_trans[vertex], &presult);
            }
        }
            break;
            
        case TRANSFORM_LOCAL_TO_TRANS: {
            for (int vertex=0; vertex < obj->num_vertices; vertex++) {
                POINT4D presult;
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex], mt, &obj->vlist_trans[vertex]);
            }
        }
            break;
        default: break;
    }
    //最后检查是否要对朝向向量进行变换;如果不进行变换,朝向向量将不再有效
    if (transform_basis) {
        //旋转物体的朝向向量
        VECTOR4D vresult;//用于储存旋转结果
        //旋转ux
        Mat_Mul_VECTOR4D_4X4(&obj->ux, mt, &vresult);
        VECTOR4D_COPY(&obj->ux, &vresult);
        //旋转uy
        Mat_Mul_VECTOR4D_4X4(&obj->uy, mt, &vresult);
        VECTOR4D_COPY(&obj->uy, &vresult);
        //旋转uz
        Mat_Mul_VECTOR4D_4X4(&obj->uz, mt, &vresult);
        VECTOR4D_COPY(&obj->uz, &vresult);
    }
}

/**************物体的局部坐标到世界坐标变换******************************************/
/*
 *  这个函数将传递进行的物体的局部/模型坐标变换为世界坐标;结果被储存在物体的变换后的顶点列表(数组vlist_trans)中
 *  遍历顶点列表,根据world_pos对所有坐标进行平移;以便将坐标变换为世界坐标;并将结果存储在数组vlist_trans[]中
 */
void
Model_To_World_OBJECT4DV1(OBJECT4DV1_PTR obj, int coord_select) {
    if (coord_select == TRANSFORM_LOCAL_TO_TRANS) {
        for (int vertex=0; vertex < obj->num_vertices; vertex++) {
            //对顶点进行平移
            VECTOR4D_ADD(&obj->vlist_local[vertex], &obj->world_pos, &obj->vlist_trans[vertex]);
        }
    }
    else {
        for (int vertex=0; vertex < obj->num_vertices; vertex++) {
            //平移顶点
            VECTOR4D_ADD(&obj->vlist_trans[vertex], &obj->world_pos, &obj->vlist_trans[vertex]);
        }
    }
}

void
Build_Model_To_World_MATRIX4X4(VECTOR4D_PTR vpos, MATRIX4X4_PTR m) {
    Mat_Init_4X4(m, 1,       0,       0,       0,
                 0,       1,       0,       0,
                 0,       0,       1,       0,
                 vpos->x, vpos->y, vpos->z, 1 );
}


/**************渲染列表的局部坐标到世界坐标变换******************************************/
/*
 大多数情况下,在对物体执行局部坐标到世界坐标变换后,才将物体的多边形插入到渲染列表中,因此,根本不需要
 对整个渲染列表执行局部坐标到世界坐标变换. 但当您不想使用OBJECT4DV1,而是将多边形直接加载到表示某些大型
 网格(如地形)的
 
 */
void
Model_To_World_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,POINT4D_PTR world_pos,int coord_select) {
    //遍历顶点列表,根据world_pos对所有坐标进行平移;以便将坐标变换为世界坐标;并将结果储存在数组vlist_trans[]中
    if (coord_select == TRANSFORM_LOCAL_TO_TRANS)  {
        for (int poly = 0; poly < rend_list->num_polys; poly++) {
            //获得当前多边形
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
            //当且仅当多边形没有被剔除和裁剪掉,同时处于活动状态且可见时才对其进行变换,但在线框引擎中,多边形是否是背面无关紧要
            if ((curr_poly==NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED ) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE) )
                continue;
            //满足条件,对其进行变换
            for (int vertex = 0; vertex < 3; vertex++)  {
                //平移顶点
                VECTOR4D_ADD(&curr_poly->vlist[vertex], world_pos, &curr_poly->tvlist[vertex]);
            }
        }
    }
    else { // TRANSFORM_TRANS_ONLY
        for (int poly = 0; poly < rend_list->num_polys; poly++){
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
            
            if ((curr_poly==NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED ) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE) )
                continue;
            
            for (int vertex = 0; vertex < 3; vertex++) {
                VECTOR4D_ADD(&curr_poly->tvlist[vertex], world_pos, &curr_poly->tvlist[vertex]);
            }
        }
    }
}



/**
 *   可以传递相机的类型,位置,朝向,裁剪面位置,视野,视口的大小,然后该函数将创建一个相机.
 *
 *  @param cam             相机对象
 *  @param cam_attr        相机属性
 *  @param cam_pos         相机的初始位置
 *  @param cam_dir         相机的初始化角度
 *  @param cam_target      UVN相机的初始化目标位置
 *  @param near_clip_z     近裁剪面
 *  @param far_clip_z      远裁剪面
 *  @param fov             视野,单位为度
 *  @param viewport_width  屏幕/视口大小
 *  @param viewport_height <#viewport_height description#>
 */
void
Init_CAM4DV1(CAM4DV1_PTR cam,int cam_attr,POINT4D_PTR cam_pos,VECTOR4D_PTR cam_dir,
             POINT4D_PTR cam_target,float near_clip_z,float far_clip_z,
             float fov,float viewport_width,float viewport_height) {
    
    cam->attr = cam_attr;              // camera attributes
    
    VECTOR4D_COPY(&cam->pos, cam_pos); // positions
    VECTOR4D_COPY(&cam->dir, cam_dir); // 欧拉相机的方向向量或角度
    // for UVN camera
    VECTOR4D_INITXYZ(&cam->u, 1,0,0);  // set to +x 轴方向
    VECTOR4D_INITXYZ(&cam->v, 0,1,0);  // set to +y
    VECTOR4D_INITXYZ(&cam->n, 0,0,1);  // set to +z
    
    if (cam_target!=NULL)
        VECTOR4D_COPY(&cam->target, cam_target); // UVN 目标位置
    else
        VECTOR4D_ZERO(&cam->target);
    
    cam->near_clip_z = near_clip_z;     // 近裁剪面
    cam->far_clip_z  = far_clip_z;      // 远裁剪面
    
    cam->viewport_width  = viewport_width;   // 视口的大小
    cam->viewport_height = viewport_height;
    
    cam->viewport_center_x = (viewport_width-1)/2; // 视口的中心
    cam->viewport_center_y = (viewport_height-1)/2;
    
    cam->aspect_ratio = (float)viewport_width/(float)viewport_height;
    
    // 将所有的矩阵都设置为单位矩阵
    MAT_IDENTITY_4X4(&cam->mcam);
    MAT_IDENTITY_4X4(&cam->mper);
    MAT_IDENTITY_4X4(&cam->mscr);
    
    // 设置相关变量
    cam->fov              = fov;
    
    // 将视平面大小设置为2x（2/ar）
    cam->viewplane_width  = 2.0;
    cam->viewplane_height = 2.0/cam->aspect_ratio;
    
    //根据fov 和视平面大小计算 视距 d
    float tan_fov_div2 = tan(DEG_TO_RAD(fov/2));
    
    cam->view_dist = (0.5)*(cam->viewplane_width)*tan_fov_div2;
    
    // 判断fov是否为90度
    if (fov == 90.0) {
        
        // 建立裁剪面!
        POINT3D pt_origin; // point on the plane
        VECTOR3D_INITXYZ(&pt_origin,0,0,0);
        
        VECTOR3D vn; // 面法线
        
        // 右裁剪面
        VECTOR3D_INITXYZ(&vn,1,0,-1); // x=z 平面
        PLANE3D_Init(&cam->rt_clip_plane, &pt_origin,  &vn, 1);
        
        // left clipping plane
        VECTOR3D_INITXYZ(&vn,-1,0,-1); // -x=z plane
        PLANE3D_Init(&cam->lt_clip_plane, &pt_origin,  &vn, 1);
        
        // top clipping plane
        VECTOR3D_INITXYZ(&vn,0,1,-1); // y=z plane
        PLANE3D_Init(&cam->tp_clip_plane, &pt_origin,  &vn, 1);
        
        // bottom clipping plane
        VECTOR3D_INITXYZ(&vn,0,-1,-1); // -y=z plane
        PLANE3D_Init(&cam->bt_clip_plane, &pt_origin,  &vn, 1);
    }
    else {                  //计算fov不为90时的裁剪面
        
        POINT3D pt_origin;  //平面上的一个点
        VECTOR3D_INITXYZ(&pt_origin,0,0,0);
        
        VECTOR3D vn;        //面法线
        /*
         由于fov不为90,因此计算起来比较复杂:首先计算表示裁剪面在平面x-z 和 y－z上的2d投影的向量
         然后计算与这两个向量垂直的向量,它就是裁剪面的法线,
         */
        
        //右裁剪面
        VECTOR3D_INITXYZ(&vn,cam->view_dist,0,-cam->viewplane_width/2.0);
        PLANE3D_Init(&cam->rt_clip_plane, &pt_origin,  &vn, 1);
        //左裁剪面,可以绕z轴反射右裁剪面的法线,来得到左裁剪面的法线
        //因为这两个裁剪面时关于z轴对称的,因此只需对x求负
        VECTOR3D_INITXYZ(&vn,-cam->view_dist,0,-cam->viewplane_width/2.0);
        PLANE3D_Init(&cam->lt_clip_plane, &pt_origin,  &vn, 1);
        //上裁剪面
        VECTOR3D_INITXYZ(&vn,0,cam->view_dist,-cam->viewplane_width/2.0);
        PLANE3D_Init(&cam->tp_clip_plane, &pt_origin,  &vn, 1);
        //下裁剪面
        VECTOR3D_INITXYZ(&vn,0,-cam->view_dist,-cam->viewplane_width/2.0);
        PLANE3D_Init(&cam->bt_clip_plane, &pt_origin,  &vn, 1);
    }
}

/**
 *  该函数根据欧拉角度计算相机变换矩阵,并存储到传入的相机对象中
 *
 *  @param cam         相机对象
 *  @param cam_rot_seq <#cam_rot_seq description#>
 */
void
Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq) {
    MATRIX4X4 mt_inv,       // 相机平移矩阵的逆矩阵
    mx_inv,                 // 相机绕x轴的旋转矩阵的逆矩阵
    my_inv,                 // 相机绕y轴的旋转矩阵的逆矩阵
    mz_inv,                 // 相机绕z轴的旋转矩阵的逆矩阵
    mrot,                   // 所有逆旋转矩阵的积
    mtmp;                   // 用于存储临时矩阵
    
    
    // step 1:根据相机位置计算相机平移矩阵的逆矩阵
    Mat_Init_4X4(&mt_inv, 1,    0,     0,     0,
                 0,    1,     0,     0,
                 0,    0,     1,     0,
                 -cam->pos.x, -cam->pos.y, -cam->pos.z, 1);
    // step 2 :创建旋转矩阵的逆矩阵
    // 要计算正规旋转矩阵的逆矩阵,可以将其转置,也可以将每个旋转角度取负
    // 首先计算3个旋转矩阵的逆矩阵
    
    //提取欧拉角度
    float theta_x = cam->dir.x;
    float theta_y = cam->dir.y;
    float theta_z = cam->dir.z;
    
    // 计算角度x的正玄 和余弦
    float cos_theta = Fast_Cos(theta_x);  // 余弦值不变 cos(-x) = cos(x)
    float sin_theta = -Fast_Sin(theta_x); // sin(-x) = -sin(x)
    
    // 创建矩阵
    Mat_Init_4X4(&mx_inv, 1,    0,         0,         0,
                 0,    cos_theta, sin_theta, 0,
                 0,   -sin_theta, cos_theta, 0,
                 0,    0,         0,         1);
    
    // 计算角度y的正弦与余弦
    cos_theta = Fast_Cos(theta_y);  // 余弦值不变 cos(-x) = cos(x)
    sin_theta = -Fast_Sin(theta_y); // sin(-x) = -sin(x)
    
    // 创建矩阵
    Mat_Init_4X4(&my_inv,cos_theta, 0, -sin_theta, 0,
                 0,         1,  0,         0,
                 sin_theta, 0,  cos_theta,  0,
                 0,         0,  0,          1);
    
    // 计算角度z的正弦 和 余弦
    cos_theta = Fast_Cos(theta_z);  // 余弦值不变 cos(-x) = cos(x)
    sin_theta = -Fast_Sin(theta_z); // sin(-x) = -sin(x)
    
    // 创建矩阵
    Mat_Init_4X4(&mz_inv, cos_theta, sin_theta, 0, 0,
                 -sin_theta, cos_theta, 0, 0,
                 0,         0,         1, 0,
                 0,         0,         0, 1);
    
    // 现在计算逆旋转矩阵的乘积
    switch(cam_rot_seq) {
        case CAM_ROT_SEQ_XYZ: {
            Mat_Mul_4X4(&mx_inv, &my_inv, &mtmp);
            Mat_Mul_4X4(&mtmp, &mz_inv, &mrot);
        } break;
            
        case CAM_ROT_SEQ_YXZ: {
            Mat_Mul_4X4(&my_inv, &mx_inv, &mtmp);
            Mat_Mul_4X4(&mtmp, &mz_inv, &mrot);
        } break;
            
        case CAM_ROT_SEQ_XZY:
        {
            Mat_Mul_4X4(&mx_inv, &mz_inv, &mtmp);
            Mat_Mul_4X4(&mtmp, &my_inv, &mrot);
        } break;
            
        case CAM_ROT_SEQ_YZX:
        {
            Mat_Mul_4X4(&my_inv, &mz_inv, &mtmp);
            Mat_Mul_4X4(&mtmp, &mx_inv, &mrot);
        } break;
            
        case CAM_ROT_SEQ_ZYX:
        {
            Mat_Mul_4X4(&mz_inv, &my_inv, &mtmp);
            Mat_Mul_4X4(&mtmp, &mx_inv, &mrot);
        } break;
            
        case CAM_ROT_SEQ_ZXY:
        {
            Mat_Mul_4X4(&mz_inv, &mx_inv, &mtmp);
            Mat_Mul_4X4(&mtmp, &my_inv, &mrot);
            
        } break;
            
        default: break;
    }
    //现在,mrot逆旋转矩阵的乘积
    //接下来将其乘以逆平移矩阵,并将结果存储到相机对象的相机变化矩阵中
    Mat_Mul_4X4(&mt_inv, &mrot, &cam->mcam);
}
/**
 *  UVN模型的相机函数如下:
 *
 *  @param cam  相机对象
 *  @param mode 参数mode指定如何计算uvn
 */
void
Build_CAM4DV1_Matrix_UVN(CAM4DV1_PTR cam, int mode) {
    MATRIX4X4 mt_inv,  // 逆相机平移矩阵
    mt_uvn,  // UVN相机变换矩阵
    mtmp;    // 用于存储临时矩阵
    
    // step 1: 根据相机位置创建逆平移矩阵
    Mat_Init_4X4(&mt_inv, 1,    0,     0,     0,
                 0,    1,     0,     0,
                 0,    0,     1,     0,
                 -cam->pos.x, -cam->pos.y, -cam->pos.z, 1);
    // step 2: 确定如何计算目标点
    if (mode == UVN_MODE_SPHERICAL) {
        // 使用球面坐标模式;需要重新计算目标点
        // 提取方位角和仰角
        float phi   = cam->dir.x; // elevation 仰角
        float theta = cam->dir.y; // heading 方位角
        
        // 计算三角函数
        float sin_phi = Fast_Sin(phi);
        float cos_phi = Fast_Cos(phi);
        
        float sin_theta = Fast_Sin(theta);
        float cos_theta = Fast_Cos(theta);
        // 计算目标点在单位球面上的位置 x,y,z
        cam->target.x = -1*sin_phi*sin_theta;
        cam->target.y =  1*cos_phi;
        cam->target.z =  1*sin_phi*cos_theta;
    }
    // 至此,有个重新计算u v n 所需的全部参数
    // Step 1: n = <目标位置 - 观察参考点>
    VECTOR4D_Build(&cam->pos, &cam->target, &cam->n);
    
    // Step 2: 将 v = <0,1,0>
    VECTOR4D_INITXYZ(&cam->v,0,1,0);
    
    // Step 3: u = (v x n)
    VECTOR4D_Cross(&cam->v,&cam->n,&cam->u);
    
    // Step 4: v = (n x u)
    VECTOR4D_Cross(&cam->n,&cam->u,&cam->v);
    
    // Step 5: 对所有向量都进行归一化
    VECTOR4D_Normalize(&cam->u);
    VECTOR4D_Normalize(&cam->v);
    VECTOR4D_Normalize(&cam->n);
    
    
    // 将uvn 代入,得到uvn旋转矩阵
    Mat_Init_4X4(&mt_uvn, cam->u.x,    cam->v.x,     cam->n.x,     0,
                 cam->u.y,    cam->v.y,     cam->n.y,     0,
                 cam->u.z,    cam->v.z,     cam->n.z,     0,
                 0,           0,            0,            1);
    
    // 将平移矩阵乘以uvn矩阵,并将结果存储到相机变换矩阵mcam中
    Mat_Mul_4X4(&mt_inv, &mt_uvn, &cam->mcam);
}

/******************************************************************************/

//1: 物体的 世界坐标 到 相机坐标 变换

/** 这个函数的速度很快,也更紧凑,它接受一个物体和一台相机作为参数,设置好物体和相机后，用该函数
 *  这是一个基于矩阵的函数,它根据传入的相机变换矩阵,将物体的世界坐标变换为相机坐标
 *  它完全不考虑多边形本身,只是对vlist_trans[]中的顶点进行变换,这是变换方法之一,您可能选择对
 *  渲染列表进行变换,因为渲染列表中的多边形表示的几何体都通过了背景剔除
 *  @param cam 相机对象
 *  @param obj 3D对象
 */
void
World_To_Camera_OBJECT4DV1(CAM4DV1_PTR cam,OBJECT4DV1_PTR obj) {
    //将物体的每个顶点变换为相机坐标;这里假设顶点已经被变换为世界坐标,且结果存储在vlist_trans[]中
    for (int vertex = 0; vertex < obj->num_vertices; vertex++) {
        //使用相机对象中的矩阵mcam对顶点进行变换
        POINT4D presult;//用于存储每次变换的结果
        //对顶点进行变换
        Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex], &cam->mcam, &presult);
        //将结果存回去
        VECTOR4D_COPY(&obj->vlist_trans[vertex], &presult);
    }
}


//2: 渲染列表的 世界坐标 到 相机坐标 变换 (执行到这个步骤,流水线已经完成了一半左右,接下来需要进行透视变换 和 屏幕 变换)
/**  基于矩阵:根据传入的相机变换矩阵将渲染列表中的每个多边形变换为相机坐标;
 *   如果在流水线的上游已经将每个物体转换为多边形并将他们插入到渲染列表中,将使用这个函数,而不是基于物体
 *   的函数对顶点进行变换,将物体转换为多边形的操作是在物体剔除,局部变换,局部坐标 到 世界坐标 变换
 *   以及北京消除之后进行的; 这样最大限度地减少了每个物体中被插入到渲染列表中的多边形数目
 *
 *    这个函数假设至少已经进行了局部坐标 到 世界坐标 变换 ;且多边形数据存储在POLYF4D1 的变换后的列表tvlist
 *
 *  执行渲染列表变换的专用函数,唯一复杂的地方是:需要遍历列表中的每个多边形;而对OBJECT4DV1进行变换
 *  的专用函数中,只需要对一个顶点列表进行变换
 *
 *  @param rend_list 渲染列表
 *  @param cam       相机对象
 */
void
World_To_Camera_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,CAM4DV1_PTR cam) {
    //就渲染列表中的每个多边形变换为相机坐标;
    for(int poly = 0; poly < rend_list->num_polys; poly++) {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        //判断多边形是够有效
        if((curr_poly == NULL) || (curr_poly->state & POLY4DV1_STATE_ACTIVE)
           || (curr_poly->state & POLY4DV1_STATE_CLIPPED)
           || (curr_poly->state &POLY4DV1_STATE_BACKFACE) ) {
            
            for (int vertex = 0 ;vertex < 3  ; vertex ++) {
                POINT4D presult;
                Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex], &cam->mcam, &presult);
                VECTOR4D_COPY(&curr_poly->tvlist[vertex], &presult);
            }
        }
    }
}

/******************************************************************************/

//物体剔除操作 culling,以避免在以后的流水线 进行变换 ()

/** 基于矩阵 根据传入的相机信息判断物体是否  在视景体内
 *  该函数接受物体的世界空间位置和平均最大半径作为参数,并根据当前的相机变换 来剔除物体(只适用于物体,因为物体变换为多边形并
 插入到主渲染列表后将不再有物体的概念)
 *
 *  @param obj        要进行剔除操作的 3d对象
 *  @param cam        剔除时使用的相机
 *  @param cull_flags 要考虑的裁剪面,其值为各种剔除标记,如果物体被剔除,将相应的设置其状态
 *
 *  @return <#return value description#>
 */
int
Cull_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam, int cull_flags) {
    //step 1: 将物体的包围球球心 变换为相机坐标
    POINT4D sphere_pos;//存储变换后的坐标
    Mat_Mul_VECTOR4D_4X4(&obj->world_pos, &cam->mcam, &sphere_pos);
    
    //step 2: 根据剔除标记对物体执行剔除操作
    if (cull_flags & CULL_OBJECT_Z_PLANE) {
        //只根据远近裁剪面来剔除物体
        //使用远近裁剪面进行测试
        if (((sphere_pos.z - obj->max_radius) > cam->far_clip_z )
            || (sphere_pos.z + obj->max_radius) < cam->near_clip_z) {
            SET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);
            return 1;
        }
    }
    
    if (cull_flags & CULL_OBJECT_X_PLANE) {
        //只根据左右裁剪面进行物体剔除,本可以使用平面方程,但使用三角形相似更容易
        //因为这是一种2D问题,如果视野为90度,问题更简单,这里假设视野不为90度
        
        //使用右裁剪面和左裁剪面 检测包围球上 最左边 和 最右边 的点
        float z_test = (0.5) * cam->viewplane_width * sphere_pos.z / cam->view_dist;
        if (((sphere_pos.x - obj->max_radius) > z_test ) //right
            || (sphere_pos.z + obj->max_radius) < z_test) { //left
            SET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);
            return 1;
        }
    }
    
    
    if (cull_flags & CULL_OBJECT_Y_PLANE) {
        //只根据上下裁剪面进行物体剔除,本可以使用平面方程,但使用三角形相似更容易
        //因为这是一种2D问题,如果视野为90度,问题更简单,这里假设视野不为90度
        
        //使用上裁剪面和下裁剪面 检测包围球上 最下边 和 最上边 的点
        float z_test = (0.5) * cam->viewplane_width * sphere_pos.z / cam->view_dist;
        if (((sphere_pos.x - obj->max_radius) > z_test ) //right
            || (sphere_pos.z + obj->max_radius) < z_test) { //left
            SET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);
            return 1;
        }
    }
    return 0;
}


//重置物体的标记
/**
 *  这个函数重置传入的状态,为变换做准备;通常是重置被剔除,被裁剪掉和背面等标记,也可以在这里做其他准备工作
 *  物体是有效的,接下来重置其各个多边形的状态标记
 *  @param obj 3D物体
 */
void
Reset_OBJECT4DV1(OBJECT4DV1_PTR obj) {
    //重置物体的被剔除标记
    RESET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);
    //重置多边形的被裁剪掉 和 背景 标记
    for (int poly = 0 ; poly < obj->num_polys; poly++) {
        POLY4DV1_PTR curr_poly = &obj->plist[poly];
        
        //判断多边形是否可见
        if (!(curr_poly->state & OBJECT4DV1_STATE_ACTIVE))
            continue; //进入下一个多边形
        RESET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
        RESET_BIT(curr_poly->state, POLY4DV1_STATE_CLIPPED);
    }
}

//1 :物体的背面消除

// 物体的背面消除（删除背向视点的多边形）;简单地说:将对物体或者渲染列表的每个多边形进行测试(将计算一个多边形到视点的方向向
// 量和多边形外向法线以及他们之间的夹角 :如果 夹角大于 90度,表示该多边形是背面,需要消除掉,是不可见的 )
// 当然,背面的概念只对单面多边形有意义,对于那些从两边都可以看到（既双面）多边形,这种测试是无意义的
// 这种测试通常是在世界空间 而不是 相机空间 进行的(因为只需要知道视点),这样可通过背面消除删除
//大量的多边形,避免对它们进行世界在坐标到相机坐标变换;这是一种不错的删除一半 几何体的方法
/**
 *  基于矩阵,根据数组vlist_trans中的顶点数据以及相机位置,消除物体的背面多边形,这里只设置多边形的背面状态
 *
 *  @param obj 3D对象
 *  @param cam 相机对象
 */
void
Remove_Backfaces_OBJECT4DV1(OBJECT4DV1_PTR obj,CAM4DV1_PTR cam) {
    //检查物体是否已被剔除
    if (obj->state & OBJECT4DV1_STATE_CULLED)
        return;
    //处理物体的每一个多边形
    for (int  poly = 0; poly < obj->num_polys; poly++) {
        POLY4DV1_PTR curr_poly = &obj->plist[poly];
        
        if (!  (curr_poly->state & POLY4DV1_STATE_ACTIVE)
            || (curr_poly->state & POLY4DV1_STATE_CLIPPED)
            || (curr_poly->state & POLY4DV1_STATE_BACKFACE)
            || (curr_poly->state & POLY4DV1_ATTR_2SIDED) )// 多边形 是双面的
            continue;
        //获取顶点列表中的顶点索引;多边形不是自包含的,而是基于物体的顶点列表
        int vindex_0 = curr_poly->vert[0];
        int vindex_1 = curr_poly->vert[1];
        int vindex_2 = curr_poly->vert[2];
        //我们将使用变换后的多边形顶点列表;因为背面消除只能在顶点被转换为世界坐标之后才能进行
        
        //需要计算多边形的面法线;顶点是按照顺时针方向排序 u=p0->p1,v=p0->p2,n = uxv
        VECTOR4D u,v,n;
        //计算u ,v
        VECTOR4D_Build(&obj->vlist_trans[vindex_0],
                       &obj->vlist_trans[vindex_1],
                       &u);
        VECTOR4D_Build(&obj->vlist_trans[vindex_0],
                       &obj->vlist_trans[vindex_1],
                       &v);
        //计算叉积
        VECTOR4D_Cross(&u, &v, &n);
        //创建指向视点的向量
        VECTOR4D  view;
        VECTOR4D_Build(&obj->vlist_trans[vindex_0], &cam->pos, &view);
        
        //计算点积
        float dp = VECTOR4D_Dot(&n, &view);
        //如果dp<0 ,则多边形是不可见的
        if (dp <= 0.0) {
            SET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
        }
    }
}



//2: 渲染列表的背面消除
/**
 *  如果在背面消除之前物体被存储在一起,或者只有游戏空间的多边形列表而没有物体,则需要对渲染列表执行
 *  背面消除.当然,如果在将组成物体的多边形插入渲染列表时执行了背面消除,则没必要对渲染列表执行 背面消除
 *  @param rend_list 渲染列表
 *  @param cam       相机对象
 */
void
Remove_Backfaces_RENDERLIST4DV1(RENDERLIST4DV1_PTR  rend_list,CAM4DV1_PTR cam) {
    for(int poly = 0;poly < rend_list->num_polys;poly++) {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        if ((curr_poly == NULL)|| ! (curr_poly->state & POLY4DV1_STATE_ACTIVE)
            || (curr_poly->state & POLY4DV1_STATE_CLIPPED)
            || (curr_poly->state & POLY4DV1_STATE_BACKFACE)
            || (curr_poly->state & POLY4DV1_ATTR_2SIDED) )// 多边形 是双面的
            continue;//move onto next poly

        //需要计算多边形的面法线;
        VECTOR4D u , v, n;
        
        VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[1], &u);
        VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[2], &v);
        
        VECTOR4D_Cross(&u, &v, &n);
        
        VECTOR4D view;
        VECTOR4D_Build(&curr_poly->tvlist[0], &cam->pos, &view);

        float dp = VECTOR4D_Dot(&n, &view);
        
        if(dp <= 0.0) SET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
    }
}



/******************************************************************************/

// 相机坐标 到 透视坐标 变换 （这阶段的变换 为: 相机坐标 -> 透视坐标）
// 简单地说: 给定点  p(X_world , Y_world , Z_world,1) 透视变换为:
// X_par = viewing_distance * X_world / Z_world
// Y_par = viewing_distance * aspect_ratio * Y_world / Z_world
// 如果将视距(d)设置为1.0,则视平面坐标将是归一化的,即坐标的范围为:
// x坐标--------  -1 到 1;  y 坐标 ----------- -1/ar 到 1/ar,其中ar为光栅化屏幕的宽高比

//1: 物体的透视变换
/**   这个函数 不是  基于 矩阵的
     这个函数根据传入的相机对象将物体的相机坐标变换为 透视坐标,它根本不关心多边形本身,而只是对vlist_trams[]中的顶点进行变换
 *  这只是执行透视变换的方法之一,你可能不采用这种方法,而是对渲染列表进行变换,因为渲染列表中的多边形表示的是透过背面消除的几何体,最后这个函数只是基于实验的目的而编写的,在真正的3d流水线中 物体不会完整地保留到这个阶段,因为物体可能只有一个多边形是可以见的,而这个函数对所有多边形都进行变换
 *  相机对象中存储了执行透视变换所需要的所有信息,这包括(视距,视平面的 长 与 宽 以及视野)
 *  waring: 这个函数中,没有检查 z值是否大于0,这种情况需要在之前的剔除 与 消除裁剪过程中进行处理,这里假设所有的顶点都是可以投影的
 *  @param obj 3D物体
 *  @param cam 相机对象
 */
void
Camera_To_Perspective_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam) {
    //将物体的每个顶点变换为 透视坐标,这里假设物体已经被变换为相机坐标,且结果存储在vlist_trans[]中
    for (int vertex = 0; vertex < obj->num_vertices; vertex++) {
        float z = obj->vlist_trans[vertex].z;
        //根据相机的观察参数对 顶点进行变换
        obj->vlist_trans[vertex].x = cam->view_dist * obj->vlist_trans[vertex].x/z;
        obj->vlist_trans[vertex].y = cam->view_dist * obj->vlist_trans[vertex].y * cam->aspect_ratio/z;
    }
}

/**
 *  矩阵来执行投影 (基于矩阵) 创建 相机坐标 到 透视坐标 变换矩阵
 *  在大多数情况下,相机的视平面是归一化的(2X2),FOV为90度,这里假设使用的是4D齐次坐标,在某个时候执行
   4D坐标到3D坐标 转换
    这种操作 可能在透视变换 之后马上进行,也可能在屏幕变换之后才进行
 *  @param cam 相机对象
 *  @param m   透视变换 矩阵
 */
void
Build_Camera_To_Perspective_MATRIX4X4(CAM4DV1_PTR cam,MATRIX4X4_PTR m) {
 Mat_Init_4X4(m,
              cam->view_dist, 0, 0, 0,
              0, cam->view_dist * cam->aspect_ratio, 0, 0,
              0, 0, 1, 1,
              0, 0, 0, 0);
}
/*  执行透视变换后,顶点的w坐标不再为1.0,必须除以w分量 将齐次坐标转换为 非齐次坐标(真正的3d坐标) */

/**
 *   这个函数将变换后的顶点列表中所有的顶点从4D齐次坐标转换为3D坐标,方法为:将分量x,y,z都除以w
 *  @param obj 3D物体
 */

void
Convert_From_Homogeneous4D_OBJECT4DV1(OBJECT4DV1_PTR obj) {
    for (int vertex = 0; vertex < obj->num_vertices; vertex++) {
        VECTOR4D_DIV_BY_W(&obj->vlist_trans[vertex]);
    }
}


//2 : 渲染列表的透视变换
/**
 *  不使用矩阵来执行 相机坐标到 透视坐标变换的函数(不基于矩阵)
 *  根据传入的相机对象,将渲染列表中的每个多边形都变换为透视坐标, 如果在流水线的上游已经将物体转换为了多边形
 *
 *  @param rend_list 渲染列表
 *  @param cam       相机对象
 */
void
Camera_To_Perspective_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam) {
    
    for (int poly =0 ; poly < rend_list->num_polys; poly++) {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        
        if ((curr_poly == NULL)
            || (curr_poly->state & POLY4DV1_STATE_ACTIVE)
            || (curr_poly->state & POLY4DV1_STATE_BACKFACE)
            || (curr_poly->state & POLY4DV1_STATE_CLIPPED)) {
            continue;
        }
        
        for (int vertex = 0; vertex < 3; vertex++) {
            float z = curr_poly->tvlist[vertex].z;
            
            curr_poly->tvlist[vertex].x = cam->view_dist * curr_poly->tvlist[vertex].x / z;
            curr_poly->tvlist[vertex].y = cam->view_dist * curr_poly->tvlist[vertex].y
            * cam->aspect_ratio / z;
        }
    }
}



/**
 *  齐次坐标转换为非齐次坐标 4d坐标 转换 为 3d坐标 ; 将分量x,y,z都除以w
 *
 *  @param rend_list 渲染列表
 */
void
Convert_From_Hogogeneous4D_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list) {
    for (int poly = 0 ; poly < rend_list->num_polys; poly++) {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        if ((curr_poly == NULL)
            || (curr_poly->state & POLY4DV1_STATE_ACTIVE)
            || (curr_poly->state & POLY4DV1_STATE_BACKFACE)
            || (curr_poly->state & POLY4DV1_STATE_CLIPPED)) {
            continue;
        }

        for (int vertex = 0 ; vertex < 3; vertex++) {
            //4d坐标 转 3d坐标
            VECTOR4D_DIV_BY_W(&curr_poly->tvlist[vertex]);
        }
    }
}



/******************************************************************************/

//透视坐标 到 屏幕(视口)坐标变换
// (从某种意义上说,这一阶段为3D流水线的最后一个阶段,它将视平面坐标进行缩放,变成屏幕坐标)
// 在这种变换中必须考虑这样一点,即大多数光栅化屏幕的原点位于左上角,y轴的方向与2d笛卡尔坐标系是相反的
//当然如果在透视变换中,视平面的大小与视口相同,则无需执行缩放操作,但大多数情况下,需要执行平移,并反转y轴
//因为在投影时,我们假设视平面的中心为 原点,其+x 轴指向右方,+y轴指向上方,而光栅屏幕的原点位于左上角
// y轴方向与此相反,因此,无论什么情况下,都需要执行某种形式的视口变换.

//1:物体的视口变换
/**
 *  不是基于 矩阵,根据传入的视口信息将物体的 透视坐标 变换 为 屏幕坐标,但完全不关心多边形本身,而只是
 *  对vlist_trans[]中的顶点进行变换,这只是执行屏幕变换的方法之一
 *  @param obj <#obj description#>
 *  @param cam <#cam description#>
 */
void
Perspective_To_Screen_OBJECT4DV1(OBJECT4DV1_PTR obj,CAM4DV1_PTR cam) {
    //将物体的每个顶点变换为屏幕坐标
    float alpha = (0.5 *cam->viewport_width - 0.5);
    float beta =  (0.5 *cam->viewport_height- 0.5);
    
    for (int vertex = 0 ; vertex < obj->num_vertices; vertex++) {
        obj->vlist_trans[vertex].x = alpha + alpha *obj->vlist_trans[vertex].x;
        obj->vlist_trans[vertex].y = beta - beta *obj->vlist_trans[vertex].y;
    }
}

//1.2:物体的视口变换 矩阵
/**
 *  这个函数创建透视坐标 到 屏幕坐标 变换矩阵
 *
 *  @param cam 相机对象
 *  @param m   矩阵对象
 */
void
Build_Perspective_To_Screen_4D_MATRIX4X4(CAM4DV1_PTR cam,MATRIX4X4_PTR m) {
    float alpha = (0.5 *cam->viewport_width - 0.5);
    float beta =  (0.5 *cam->viewport_height- 0.5);
    Mat_Init_4X4(m, alpha,   0,     0,    0,
                 0,      -beta,  0,    0,
                 alpha,   beta,  1,    0,
                 0,       0,     0,    1);

}


//2 :渲染列表 到 视口变换
/**
 *  非矩阵;这个函数只对 坐标进行缩放,并反转y轴
 *
 *  @param rend_list 渲染列表
 *  @param cam       相机对象
 */
void
Perspective_To_Screen_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list ,CAM4DV1_PTR cam) {
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        if ((curr_poly==NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED ) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE) )
            continue; // move onto next poly
        
        float alpha = (0.5*cam->viewport_width-0.5);
        float beta  = (0.5*cam->viewport_height-0.5);
        
        for (int vertex = 0; vertex < 3; vertex++){
            curr_poly->tvlist[vertex].x = alpha + alpha*curr_poly->tvlist[vertex].x;
            curr_poly->tvlist[vertex].y = beta  - beta *curr_poly->tvlist[vertex].y;
        }
    }
}



/******************************************************************************/

// 合并透视变换 和  屏幕变换 （非矩阵）
//一般而言,执行相机坐标 到 屏幕坐标 变换时,最快捷的方式是创建一个手工完成这项工作并将x和y坐标 除以z的函数

//1 : 物体的 相机坐标 到 屏幕坐标 变换
/**
 *  这个函数减少了更多地数学计算,我们无需根据视口的大小相应地缩放视平面坐标,因为物体被投影到大小为
 *  viewport_width * viewport_height的视平面上,这进一步减少了计算量
 *
 *  @param obj <#obj description#>
 *  @param cam <#cam description#>
 */
void
Camera_To_Perspective_Screen_OBJECT4DV1(OBJECT4DV1_PTR obj,CAM4DV1_PTR cam) {
    float alpha = (0.5*cam->viewport_width-0.5);
    float beta  = (0.5*cam->viewport_height-0.5);
    
    for (int vertex = 0; vertex < obj->num_vertices; vertex++) {
    
        float z = obj->vlist_trans[vertex].z;
        obj->vlist_trans[vertex].x = cam->view_dist*obj->vlist_trans[vertex].x/z;
        obj->vlist_trans[vertex].y = cam->view_dist*obj->vlist_trans[vertex].y/z;
        obj->vlist_trans[vertex].x =  obj->vlist_trans[vertex].x + alpha;
        obj->vlist_trans[vertex].y = -obj->vlist_trans[vertex].y + beta;
        
    }
}

//2 渲染列表的 相机坐标 到 屏幕坐标 变换  （到这个阶段,3D流水线的大部分工作已经完成,我们可以一次性执行透视变换和屏幕变换,以节省时间和计算量,大多数情况下，我们采用这种合二为一的方法,至少使用手工方法来执行透视变换和屏幕变换）
//  将坐标变换为3D的，因为不想将4D坐标到3D坐标转换的工作 留到流水线的后续阶段进行
/**
 *  对渲染列表一次执行 相机坐标到 屏幕坐标 变换 （不是基于矩阵的）
 *
 *  @param rend_list <#rend_list description#>
 *  @param cam       <#cam description#>
 */
void
Camera_To_Perspective_Screen_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,CAM4DV1_PTR cam) {
    for (int poly = 0; poly < rend_list->num_polys; poly++) {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        
        if ((curr_poly==NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED ) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE) )
            continue; // move onto next poly
        
        float alpha = (0.5*cam->viewport_width-0.5);
        float beta  = (0.5*cam->viewport_height-0.5);
        
        for (int vertex = 0; vertex < 3; vertex++) {
        
            float z = curr_poly->tvlist[vertex].z;
            
            curr_poly->tvlist[vertex].x = cam->view_dist*curr_poly->tvlist[vertex].x/z;
            curr_poly->tvlist[vertex].y = cam->view_dist*curr_poly->tvlist[vertex].y/z;
            curr_poly->tvlist[vertex].x =  curr_poly->tvlist[vertex].x + alpha;
            curr_poly->tvlist[vertex].y = -curr_poly->tvlist[vertex].y + beta;
        }
    }
}


/**
 *  渲染3D世界 (现在的3D流水线)
  阶段0 : 加载并放置物体
  阶段1 : 局部物体坐标 到 世界坐标 变换
  阶段2 : 物体剔除  和 背面消除
  阶段3 : 世界坐标 到 相机坐标 变换
  阶段4 : 相机坐标 到 透视坐标 变换
  阶段5 : 透视坐标 到 屏幕坐标 变换
  阶段6 : 渲染几何体(线框引擎)
 */















