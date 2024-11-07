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

// Ԥ����ָ��
// ��ֹͷ�ļ�����ΰ�����Ҳ��Ϊ��������ʿ����ͷ�ļ���ʿ������
#ifndef GLOBALMATCH_COMMON_H
#define GLOBALMATCH_COMMON_H

#include <pcl/common/common.h>

// typedef�ؼ�������Ϊ�������Ͷ���һ���µ����ƣ�������

typedef pcl::PointXY Point2D;
typedef pcl::PointXYZ Point3D;
typedef pcl::PointCloud <Point2D> Cloud2D;
typedef pcl::PointCloud <Point3D> Cloud3D;

// std::vector������C++��׼ģ��⣨STL���е�һ�������������ܹ��洢�ɱ�������ͬ����Ԫ�ء�
// <Cloud2D, Eigen::aligned_allocator<Cloud2D>>������std::vectorģ��Ĳ������֡�
// Cloud2DΪ���Ͳ���
// Eigen::aligned_allocator<Cloud2D>Ϊ��ѡ�ķ��������Ͳ���,��Eigen�����ṩ��һ����������ȷ��������ڴ��Ƕ���ġ�

// ���Cloud2D��һ������Eigen�����ݽṹ�����ͣ����磬��������Eigen�ľ������������
// ���������㽫���洢��STL��������std::vector���У���ôʹ��Eigen::aligned_allocator<Cloud2D>ͨ���Ǳ�Ҫ�ġ�
// ������ΪEigen�����ݽṹҪ���ض����ڴ���룬���Ż����ܲ���������ʱ����
typedef std::vector <Cloud2D, Eigen::aligned_allocator<Cloud2D>> Cloud2DVector;
typedef std::vector <Cloud3D, Eigen::aligned_allocator<Cloud3D>> Cloud3DVector;
typedef pcl::PointCloud <pcl::Normal> CloudNormal;

#endif //GLOBALMATCH_COMMON_H
