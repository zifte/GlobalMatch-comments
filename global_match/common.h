/*
 * Software License Agreement (Apache License)
 *
 *  Copyright (C) 2023, Xufei Wang (tjwangxufei@tongji.edu.cn),
 *                      Zexin Yang (zexinyang@tongji.edu.cn),
 *                      Liangliang Nan (liangliang.nan@gmail.com).
 *  All rights reserved.
 *
 *  This file is part of GlobalMatch (https://github.com/zexinyang/GlobalMatch),
 *  which implements the point cloud registration method described in the following paper:
 *  -----------------------------------------------------------------------------------------------------------
 *  GlobalMatch: Registration of forest terrestrial point clouds by global matching of relative stem positions.
 *  Xufei Wang, Zexin Yang, Xiaojun Cheng, Jantien Stoter, Wenbing Xu, Zhenlun Wu, and Liangliang Nan.
 *  ISPRS Journal of Photogrammetry and Remote Sensing. Vol. 197, 71-86, 2023.
 *  -----------------------------------------------------------------------------------------------------------
 *  We kindly ask you to cite the above paper if you use (part of) the code or ideas in your academic work.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */

// 预处理指令
// 防止头文件被多次包含（也称为“包含卫士”或“头文件卫士”）。
#ifndef GLOBALMATCH_COMMON_H
#define GLOBALMATCH_COMMON_H

#include <pcl/common/common.h>

// typedef关键字用于为数据类型定义一个新的名称（别名）

typedef pcl::PointXY Point2D;
typedef pcl::PointXYZ Point3D;
typedef pcl::PointCloud <Point2D> Cloud2D;
typedef pcl::PointCloud <Point3D> Cloud3D;

// std::vector：这是C++标准模板库（STL）中的一个序列容器，能够存储可变数量的同类型元素。
// <Cloud2D, Eigen::aligned_allocator<Cloud2D>>：这是std::vector模板的参数部分。
// Cloud2D为类型参数
// Eigen::aligned_allocator<Cloud2D>为可选的分配器类型参数,是Eigen库中提供的一个分配器，确保分配的内存是对齐的。

// 如果Cloud2D是一个包含Eigen库数据结构的类型（例如，它包含了Eigen的矩阵或向量），
// 并且您打算将它存储在STL容器（如std::vector）中，那么使用Eigen::aligned_allocator<Cloud2D>通常是必要的。
// 这是因为Eigen的数据结构要求特定的内存对齐，以优化性能并避免运行时错误。
typedef std::vector <Cloud2D, Eigen::aligned_allocator<Cloud2D>> Cloud2DVector;
typedef std::vector <Cloud3D, Eigen::aligned_allocator<Cloud3D>> Cloud3DVector;
typedef pcl::PointCloud <pcl::Normal> CloudNormal;

#endif //GLOBALMATCH_COMMON_H
