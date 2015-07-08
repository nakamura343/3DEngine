# 3DMathLib

(2015.07.07 新增) :开始创建3D引擎的每个步骤--从数据表示到最后的光栅化,不教怎么使用OpenGL图形API,只是检验自己是否真正的理解了3D图形学,必须知道如何手动编写类似OpenGL,这是我追求的目标:在任何计算机上编写出基于软件的3D图形引擎,在具备这种能力后,学习3D API将是小菜一碟.



当你编写3d引擎的代码时,需要的一些函数和宏很可能我已经编写好了。

在编写软件图形引擎的时候,唯一的问题就是表示方法,即应应用3D点还是齐次4D点?我们拭目以待,但至少知道如何处理这两者.

实际上,我将使用4x3和4x4矩阵,所以从技术上来说,我将使用4D齐次点,并假设w=1或将w作为第4个分量



3DMathLib  C/C++ 因为我要实现一套  3d软件渲染APi 所以重写自己的3d数学库

主要处理点，向量，直线，矩阵，四元数,参数化直线,3D平面,极坐标,柱面坐标,球面坐标和定点数

命名规则对于每个人来说都是问题，因为您希望函数有个合理的名称，但又不想输入太多的字符。
在大部分情况下，我使用类/结构名作为函数名称的一部分,加上函数功能的英文描述,并使用下划线
来连接。

eg: void Mat_Mul_VECTOR3D_3X3(VECTOR3D_PTR va,
                              MATRIX3X3_PTR mb,
                              VECTOR3D_PTR vprod);

上述函数声明非常明确的说明了函数的功能，OpenGL采用了类似的命名规则

建立数学库

1:数学常量

2:数据结构

3:宏和内联函数

4:函数原型

5:全局变量

6:坐标系支持

7:向量支持

8:向量支持

9:矩阵支持

10:2d和3d参数化直线支持

11:3d平面支持

12:四元数支持

13:定点数支持(纯学习)

14:浮点数支持
