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

#ifndef GLOBALMATCH_STEM_MATCHING_H
#define GLOBALMATCH_STEM_MATCHING_H

#include "common.h"
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/distances.h>
#include <pcl/PolygonMesh.h>
#include <pcl/registration/transformation_estimation_2D.h>
#include <algorithm> // 这个头文件包含了各种算法的实现，比如排序（std::sort）、查找（std::find）、变换（std::transform）、填充（std::fill）等。这些算法通常用于操作容器中的元素。
#include <vector> // std::vector是一个动态数组，能够在运行时动态地增加或减少元素。它提供了随机访问迭代器，允许在常数时间内访问任何元素。
#include <set> // std::set是一个有序集合，它自动对其元素进行排序，并且不允许重复元素。std::set基于红黑树实现，提供了快速的查找、插入和删除操作。
#include <fstream> // 这个头文件包含了文件输入/输出（I/O）的类，比如std::ifstream（用于从文件读取数据）和std::ofstream（用于向文件写入数据）。这些类提供了对文件的基本读写操作。
#include <numeric> // 这个头文件包含了一些数值算法，比如部分和（std::accumulate）、相邻差（std::adjacent_difference）等。这些算法通常用于对容器中的元素进行数值计算。

class Matching {
public:
    Matching() : knn_(20),
                 max_stems_for_exhaustive_search_(50),
                 edge_diff_(0.05),
                 stem_positions_src_(new Cloud3D),
                 stem_positions_tgt_(new Cloud3D) {
    }

    inline void
    setPairwiseStemPositions(const Cloud3D::ConstPtr& stem_positions_src,
                             const Cloud3D::ConstPtr& stem_positions_tgt) {
        stem_positions_src_ = stem_positions_src;
        stem_positions_tgt_ = stem_positions_tgt;
    }

    void
    estimateTransformation(Eigen::Matrix4f& transform);

    // size_t 是一个无符号整数类型，通常用于表示对象的大小或数组中的元素数量。
    // 在这个上下文中，它用于表示匹配的数量。
    inline size_t 
    getNumberOfMatches() {
        return stem_matches_.size(); // 顶层容器的元素数量
    };

    // TODO: add setters and getters

private:
    struct VertexSide {
        int vertex; // 点的索引
        float side;
    };
    typedef std::vector<VertexSide> Triangle;

    // 成员函数后面的const关键字表示该成员函数是一个常量成员函数（const member function）
    // 函数末尾的 const 表示该函数不会修改当前对象的任何成员变量。
    void
    constructTriangles(const Cloud3D::ConstPtr& stem_positions,
                       std::vector<Triangle>& triangles) const;

    bool
    satisfyLocalConsistency(const Point3D& feature_src,
                            const Point3D& feature_tgt) const;

    bool
    satisfyGlobalConsistency(const std::vector<int>& pair_initial,
                             const std::vector<int>& pair_candidate);

    void
    localMatching();

    void
    globalMatching();

    /**
     * @brief Accelerate global matching by growing only a limited number of (e.g., 10k)
     * randomly sampled groups of triangle pairs. Use this function if your point cloud contains
     * a large number of trees (e.g., more than 10k trees per point cloud).
     */
    void
    randomGlobalMatching();

    Cloud3D::ConstPtr stem_positions_src_, stem_positions_tgt_;
    std::vector<Triangle> triangles_src_, triangles_tgt_;
    std::vector<std::pair<int, int>> locally_matched_pairs_;
    std::vector<int> globally_matched_pairs_;
    pcl::Correspondences stem_matches_;

    int knn_;
    int max_stems_for_exhaustive_search_;
    float edge_diff_;
};


#endif //GLOBALMATCH_STEM_MATCHING_H
