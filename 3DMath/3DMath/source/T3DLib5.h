//
//  T3DLib5.h
//  3DMath
//
//  Created by TTc on 15/4/10.
//  Copyright (c) 2015年 TTc. All rights reserved.
//

#ifndef ___DMath__T3DLib5__
#define ___DMath__T3DLib5__

#include <stdio.h>
#include "T3DLib4.h"
// TYPES //////////////////////////////////////////////////////

// a polygon based on an external vertex list
// 基于顶点列表的多边形
typedef struct POLY4DV1_TYP{
    int state;  // 状态信息
    int attr;   //多边形的物理属性
    int color;  //多边形的颜色
    
    POINT4D_PTR vlist;//顶点列表  the vertex list itself
    int vert[3];      //顶点列表中元素的索引 the indices into the vertex list
} POLY4DV1, *POLY4DV1_PTR;

//自包含的多边形数据结构,供渲染列表使用
// a self contained polygon used for render list
typedef struct POLYF4DV1_TYP{
    int state;  //state information
    int attr;   // physical attributes of polygon
    int color;  //color of polygon
    
    POINT4D vlist[3];   //该三角形的顶点
    POINT4D tvlist[3];  //交换后的顶点
    
    POLYF4DV1_TYP *next; //指向渲染列表中下一个多边形的指针
    POLYF4DV1_TYP *prev; //指向渲染列表中前一个多边形的指针
    
} POLYF4DV1, *POLYF4DV1_PTR;

//多边形盒多边形面的属性
#define POLY4DV1_ATTR_2SIDED              0x0001
#define POLY4DV1_ATTR_TRANSPARENT         0x0002
#define POLY4DV1_ATTR_8BITCOLOR           0x0004
#define POLY4DV1_ATTR_RGB16               0x0008
#define POLY4DV1_ATTR_RGB24               0x0010

#define POLY4DV1_ATTR_SHADE_MODE_PURE       0x0020
#define POLY4DV1_ATTR_SHADE_MODE_CONSTANT   0x0020 // (alias)
#define POLY4DV1_ATTR_SHADE_MODE_FLAT       0x0040
#define POLY4DV1_ATTR_SHADE_MODE_GOURAUD    0x0080
#define POLY4DV1_ATTR_SHADE_MODE_PHONG      0x0100
#define POLY4DV1_ATTR_SHADE_MODE_FASTPHONG  0x0100 // (alias)
#define POLY4DV1_ATTR_SHADE_MODE_TEXTURE    0x0200

// states of polygons and faces
// 多边形和面的状态值
#define POLY4DV1_STATE_ACTIVE             0x0001
#define POLY4DV1_STATE_CLIPPED            0x0002
#define POLY4DV1_STATE_BACKFACE           0x0004

// defines for objects version 1
#define OBJECT4DV1_MAX_VERTICES           1024  // 64
#define OBJECT4DV1_MAX_POLYS              1024 // 128

//states for objects
#define OBJECT4DV1_STATE_ACTIVE           0x0001
#define OBJECT4DV1_STATE_VISIBLE          0x0002
#define OBJECT4DV1_STATE_CULLED           0x0004


// render list defines
#define RENDERLIST4DV1_MAX_POLYS          32768// 16384


#define TRANSFORM_LOCAL_ONLY       0  //对局部/模型顶点列表进行变换

#define TRANSFORM_TRANS_ONLY       1   //对变换后的顶点列表进行变换

#define TRANSFORM_LOCAL_TO_TRANS   2   //对局部顶点列表进行变换,并将结果存储在变换后的顶点列表中

//general culling flags 剔除标记状态
#define CULL_OBJECT_X_PLANE           0x0001 // 根据左右裁剪面进行剔除
#define CULL_OBJECT_Y_PLANE           0x0002 // cull on the y clipping planes
#define CULL_OBJECT_Z_PLANE           0x0004 // cull on the z clipping planes
#define CULL_OBJECT_XYZ_PLANES        (CULL_OBJECT_X_PLANE | CULL_OBJECT_Y_PLANE | CULL_OBJECT_Z_PLANE)


// defines for camera rotation sequences
#define CAM_ROT_SEQ_XYZ  0
#define CAM_ROT_SEQ_YXZ  1
#define CAM_ROT_SEQ_XZY  2
#define CAM_ROT_SEQ_YZX  3
#define CAM_ROT_SEQ_ZYX  4
#define CAM_ROT_SEQ_ZXY  5

// 定义相机旋转顺序对应的值
#define CAM_PROJ_NORMALIZED        0x0001
#define CAM_PROJ_SCREEN            0x0002
#define CAM_PROJ_FOV90             0x0004

#define CAM_MODEL_EULER            0x0008
#define CAM_MODEL_UVN              0x0010

#define UVN_MODE_SIMPLE            0    //低级简单模型,使用目标位置和观察参考点
#define UVN_MODE_SPHERICAL         1    //球面坐标模式,分量x和y被用作观察向量的方位角和仰角


// an object based on a vertex list and list of polygons
typedef struct OBJECT4DV1_TYP {

    int  id;           // 物体的数字ID
    char name[64];     // 物体的字符名称
    int  state;        // 物体的状态
    int  attr;         // 物体的属性
    float avg_radius;  // 物体的平均半径,用于碰撞检测
    float max_radius;  // 物体的最大半径
    
    POINT4D world_pos;  // 物体在世界坐标系中的位置
    
    VECTOR4D dir;       // 物体在局部坐标系中的旋转角度;用户定义的坐标系或单位方向向量
    
    VECTOR4D ux,uy,uz;  // 局部坐标轴,用于储存物体的朝向;旋转期间将被自动更新
    
    int num_vertices;   // 物体的顶点数
    
    POINT4D vlist_local[OBJECT4DV1_MAX_VERTICES]; // 用于储存顶点局部坐标的数组
    POINT4D vlist_trans[OBJECT4DV1_MAX_VERTICES]; // 储存变换后的顶点坐标的数组
    
    int num_polys;        // 物体网格的多边形数
    POLY4DV1 plist[OBJECT4DV1_MAX_POLYS];  // 多边形数组
    
} OBJECT4DV1, *OBJECT4DV1_PTR;


//储存渲染列表的对象,这样可以同时有多个渲染列表
typedef struct RENDERLIST4DV1_TYP {
    int state;           //渲染列表的状态
    int attr;            //渲染列表的属性
    
    //渲染列表是一个指针数组;其中每个指针指向一个自包含的,可渲染的多边形面
    POLYF4DV1_PTR poly_ptrs[RENDERLIST4DV1_MAX_POLYS];
    //为避免每帧都为多边形分配和释放存储空间
    POLYF4DV1 poly_data[RENDERLIST4DV1_MAX_POLYS];

    int num_polys;      //渲染列表中包含的多边形数目

} RENDERLIST4DV1, *RENDERLIST4DV1_PTR;


typedef struct CAM4DV1_TYP {
    int state;      //相机状态
    int attr ;      //相机属性
    
    POINT4D pos;    //相机在世界坐标系中的位置
    POINT4D dir;    //欧拉角度或UVN相机模型的注视方向
    
    VECTOR4D u;     //uvn 相机模型的朝向向量
    VECTOR4D v;
    VECTOR4D n;
    
    VECTOR4D target;    //uvn模型的目标位置
    
    float view_dist;    //水平视距和垂直视距,在透视变换中需要使用它们
    
    float fov;      //水平方向和垂直方向的视野
    
    //3d裁剪面
    //如果视野不是90度,3d裁剪面方程将为一般性平面方程
    float near_clip_z;  // 近裁剪面
    float far_clip_z;   // 远裁剪面
    
    PLANE3D rt_clip_plane;  // 右裁剪面
    PLANE3D lt_clip_plane;  // 左裁剪面
    PLANE3D tp_clip_plane;  // 上裁剪面
    PLANE3D bt_clip_plane;  // 下裁剪面
    
    float viewplane_width;     // 视平面的宽度和高度
    float viewplane_height;    // 对于归一化投影,为2x2;否则大小与视口或屏幕窗口相同
     // 屏幕和视口是同义词
    float viewport_width;      //屏幕/视口的大小
    float viewport_height;
    float viewport_center_x;   //视口的中心
    float viewport_center_y;
    
    // 宽高比
    float aspect_ratio;        //屏幕的宽高比
 
    MATRIX4X4 mcam;          //用于存储世界坐标到相机坐标变换矩阵
    MATRIX4X4 mper;          //用于存储相机坐标到透视坐标变换矩阵
    MATRIX4X4 mscr;          //用于存储透视坐标到屏幕坐标变换矩阵
    
} CAM4DV1, *CAM4DV1_PTR;

//functions

void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list);

void Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,MATRIX4X4_PTR mt,int coord_select);
//应用变换矩阵,对世界坐标进行变换
void Transform_OBJECT4DV1(OBJECT4DV1_PTR obj,MATRIX4X4_PTR mt, int coord_select, int transform_basis);

void Model_To_World_OBJECT4DV1(OBJECT4DV1_PTR obj, int coord_select) ;

void Build_Model_To_World_MATRIX4X4(VECTOR4D_PTR vpos, MATRIX4X4_PTR m) ;

void Model_To_World_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,POINT4D_PTR world_pos,int coord_select);


void Init_CAM4DV1(CAM4DV1_PTR cam, int attr, POINT4D_PTR cam_pos,
                  VECTOR4D_PTR cam_dir, VECTOR4D_PTR cam_target,
                  float near_clip_z, float far_clip_z, float fov,
                  float viewport_width,  float viewport_height);
//生成欧拉相机 变换矩阵
void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq);
//UVN相机
void Build_CAM4DV1_Matrix_UVN(CAM4DV1_PTR cam, int mode);
//物体世界坐标 到 相机坐标变换
void World_To_Camera_OBJECT4DV1(CAM4DV1_PTR cam,OBJECT4DV1_PTR obj);
//渲染列表的 世界坐标 到 相机坐标 变换
void World_To_Camera_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,CAM4DV1_PTR cam);


//物体剔除操作 culling,以避免在以后的流水线 进行变换
int Cull_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam, int cull_flags) ;
#endif /* defined(___DMath__T3DLib5__) */
