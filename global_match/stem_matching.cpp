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
#include "stem_matching.h"

#include <random>
// ICRA 2015 octree
#include "octree_unibn.hpp"

void
Matching::constructTriangles(const Cloud3D::ConstPtr& stem_positions,
                             std::vector<Triangle>& triangles) const {
    const auto& num_stem = stem_positions->size(); // & 可以删除
    // construct triangle vertices
    // std::set 保证集合中的每个 std::vector<int> 都是唯一的。
    // 利用 std::set 去除重复组合，确保每个三角形顶点组合唯一
    std::set<std::vector<int>> set_triangle_vertices;
    if (num_stem <= max_stems_for_exhaustive_search_) { // exhaustive search
        for (int i_stem = 0; i_stem < num_stem - 2; ++i_stem) {
            for (int j_stem = i_stem + 1; j_stem < num_stem - 1; ++j_stem) {
                for (int k_stem = j_stem + 1; k_stem < num_stem; ++k_stem) {
                    const auto& vertex_a = i_stem;
                    const auto& vertex_b = j_stem;
                    const auto& vertex_c = k_stem;
                    std::vector<int> triangle({vertex_a, vertex_b, vertex_c});
                    // 将三角形的顶点索引排序，以确保相同顶点的组合顺序一致，从而实现去重效果。
                    sort(triangle.begin(), triangle.end());
                    set_triangle_vertices.insert(triangle);
                }
            }
        }
    } else {  // neighbor-based search
        // search neighbors
        std::vector<std::vector<int>> neighbors;
        pcl::KdTreeFLANN<Point3D> kdtree;
        kdtree.setInputCloud(stem_positions);
        for (int i = 0; i < num_stem; ++i) {
            std::vector<int> ptid_nearest(knn_ + 1); // containing point itself
            std::vector<float> squared_dist_nearest(knn_ + 1);
            if (kdtree.nearestKSearch(stem_positions->points[i],
                                      knn_ + 1,
                                      ptid_nearest,
                                      squared_dist_nearest) > 0)
                neighbors.push_back(ptid_nearest);
        }
        // construct triangle vertices
        for (int i_stem = 0; i_stem < num_stem; ++i_stem) {
            // j_stem starts from 1 here
            // because the nearest neighbor of a point is the point itself
            for (int j_stem = 1; j_stem < knn_; ++j_stem) {
                for (int k_stem = j_stem + 1; k_stem < knn_ + 1; ++k_stem) {
                    const auto& vertex_a = i_stem;
                    const auto& vertex_b = neighbors[i_stem][j_stem];
                    const auto& vertex_c = neighbors[i_stem][k_stem];
                    std::vector<int> triangle({vertex_a, vertex_b, vertex_c});
                    // 将三角形的顶点索引排序，以确保相同顶点的组合顺序一致，从而实现去重效果。
                    sort(triangle.begin(), triangle.end());
                    set_triangle_vertices.insert(triangle);
                }
            }
        }
    }

    auto construct_a_triangle = [](const Cloud3D::ConstPtr& stem_positions,
                                   const int vertex_a,
                                   const int vertex_b,
                                   const int vertex_c) {
        float side_ab = pcl::euclideanDistance(
                stem_positions->points[vertex_a],
                stem_positions->points[vertex_b]);
        float side_ac = pcl::euclideanDistance(
                stem_positions->points[vertex_a],
                stem_positions->points[vertex_c]);
        float side_bc = pcl::euclideanDistance(
                stem_positions->points[vertex_b],
                stem_positions->points[vertex_c]);
        Triangle triangle{  // typedef std::vector<VertexSide> Triangle;
                {vertex_a, side_bc},
                {vertex_b, side_ac},
                {vertex_c, side_ab}
        };
        return triangle;
    };

    // construct triangle sides
    const auto& num_triangles = set_triangle_vertices.size();
    std::vector<std::vector<int>> vec_triangle_vertices(num_triangles);
    // 将 set_triangle_vertices 中的所有元素赋值给 vec_triangle_vertices，
    // 并且会按照 set_triangle_vertices 中的顺序（通常是升序，因为 std::set 自动对元素进行排序）
    // 插入到 vec_triangle_vertices 中。
    vec_triangle_vertices.assign(set_triangle_vertices.begin(), set_triangle_vertices.end());
    triangles.resize(num_triangles); // output
    for (int i = 0; i < num_triangles; ++i) {
        const auto& vertex_a = vec_triangle_vertices[i][0];
        const auto& vertex_b = vec_triangle_vertices[i][1];
        const auto& vertex_c = vec_triangle_vertices[i][2];
        triangles[i] = construct_a_triangle(stem_positions, vertex_a, vertex_b, vertex_c);
    }

    // sort vertices and sides by the length of sides in descending order
    for (auto& triangle: triangles) {
        // no need to worry about isosceles triangles
        /*bool sortBySide(const VertexSide & a, const VertexSide & b) {
            return a.side > b.side;
        }*/

        // lambda 表达式
        // 在 lambda 表达式中，[] 表示捕获列表。捕获列表用于指定 lambda 表达式可以访问的
        // 外部变量。它可以为空（[]），表示不捕获任何变量；
        // 也可以包含变量名（按值捕获）或 & 变量名（按引用捕获），
        // 或者使用 = 或 & 来捕获所有外部变量（分别按值或按引用）。
        // 
        // std::sort 会把 triangle.begin() 和 triangle.end() 之间的每一对元素传递给
        // compareBySide，并将这两个元素作为参数传递给 compareBySide 函数。
        auto sort_by_side = [](const VertexSide& a, const VertexSide& b) { return a.side > b.side; };
        // 按照边的长度从大到小排序
        // std::sort 在排序过程中比较的是相邻的元素，
        // 但通过递归调用和元素交换，最终会将所有元素按照要求的顺序排列。
        std::sort(triangle.begin(), triangle.end(), sort_by_side);
    }

    // sort vertices and sides in clockwise order starting from the vertex opposite to the longest side
    for (auto& triangle: triangles) {
        const auto& ver_a = stem_positions->points[triangle[0].vertex];
        const auto& ver_b = stem_positions->points[triangle[1].vertex];
        const auto& ver_c = stem_positions->points[triangle[2].vertex];
        // AB \cross_product AC, negative: clockwise, positive: counterclockwise
        float cp = (ver_b.x - ver_a.x) * (ver_c.y - ver_a.y) - (ver_b.y - ver_a.y) * (ver_c.x - ver_a.x);
        if (cp > 0) std::swap(triangle[1], triangle[2]);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
Matching::satisfyLocalConsistency(const Point3D& feature_src,
                                  const Point3D& feature_tgt) const {
    return abs(feature_src.x - feature_tgt.x) <= edge_diff_ &&
           abs(feature_src.y - feature_tgt.y) <= edge_diff_ &&
           abs(feature_src.z - feature_tgt.z) <= edge_diff_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
Matching::satisfyGlobalConsistency(const std::vector<int>& pair_initial,
                                   const std::vector<int>& pair_candidate) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            float dist_src = pcl::euclideanDistance(
                    stem_positions_src_->points[pair_initial[i]],
                    stem_positions_src_->points[pair_candidate[j]]);
            float dist_tgt = pcl::euclideanDistance(
                    stem_positions_tgt_->points[pair_initial[i + 3]],
                    stem_positions_tgt_->points[pair_candidate[j + 3]]);

            if (abs(dist_tgt - dist_src) > edge_diff_)
                return false;
        }
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Matching::localMatching() {
    const auto& num_tri_src = triangles_src_.size();
    const auto& num_tri_tgt = triangles_tgt_.size();
    std::vector<int> tgt_to_src(num_tri_tgt);

    Cloud3D::Ptr features_src(new Cloud3D), features_tgt(new Cloud3D);
    features_src->resize(num_tri_src);
    for (size_t i_src = 0; i_src < num_tri_src; ++i_src) {
        features_src->points[i_src].x = triangles_src_[i_src][0].side;
        features_src->points[i_src].y = triangles_src_[i_src][1].side;
        features_src->points[i_src].z = triangles_src_[i_src][2].side;
    }
    features_tgt->resize(num_tri_tgt);
    for (size_t i_tgt = 0; i_tgt < num_tri_tgt; ++i_tgt) {
        features_tgt->points[i_tgt].x = triangles_tgt_[i_tgt][0].side;
        features_tgt->points[i_tgt].y = triangles_tgt_[i_tgt][1].side;
        features_tgt->points[i_tgt].z = triangles_tgt_[i_tgt][2].side;
    }

    pcl::KdTreeFLANN<Point3D> kdtree;
    kdtree.setInputCloud(features_src);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int i_tgt = 0; i_tgt < num_tri_tgt; ++i_tgt) {
        std::vector<int> ptid_nearest(1); // 大小为 1 
        std::vector<float> squared_dist_nearest(1); // 大小为 1 
        kdtree.nearestKSearch(features_tgt->points[i_tgt],
                              1,
                              ptid_nearest,
                              squared_dist_nearest);
        if (satisfyLocalConsistency(
                features_src->points[ptid_nearest[0]],
                features_tgt->points[i_tgt]))
            tgt_to_src[i_tgt] = ptid_nearest[0];
        else
            tgt_to_src[i_tgt] = -1;
    }

    for (int i = 0; i < num_tri_tgt; ++i) {
        if (tgt_to_src[i] != -1)
            // emplace_back 是 std::vector 的一个成员函数，用于在向量末尾插入一个新元素。
            // 在这里，它将 tgt_to_src[i]（源三角形的索引）和 i（目标三角形的索引）作为一个配对插入到
            // locally_matched_pairs_ 中。
            locally_matched_pairs_.emplace_back(tgt_to_src[i], i); // source to target
        // emplace_back和push_back的区别:
        // 如果你已经有一个对象，并且想要将它加入容器，使用 push_back。
        // 如果你只想通过构造函数在容器末尾创建一个新对象，则应该使用 emplace_back。
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Matching::globalMatching() {
    const auto& num_pairs = locally_matched_pairs_.size();
    std::vector<std::vector<int>> groups(num_pairs); // 二维向量
    std::vector<int> group_size(num_pairs); // 一维向量
#if defined(_OPENMP)
    // 这是一个OpenMP的编译器指令（pragma），用于指示接下来的for循环应该被并行化。
    // parallel关键字告诉编译器这个循环应该并行执行，
    // for则指定了接下来的循环结构将被并行化。
    // schedule(static)指定了迭代被静态地分配给各个线程，每个线程处理固定数量的迭代。
#pragma omp parallel for schedule(static)
#endif
    for (int i_grp = 0; i_grp < num_pairs; ++i_grp) {// 每一对都作为初始组，使用剩下的对计算全局偏差，之后保留符合条件的最多组数
        groups[i_grp].push_back(i_grp);
        // initial triangle pair
        std::vector<int> pair_init(6);
        for (int i = 0; i < 3; ++i) // source triangle
            // 这里用.first访问了第一个元素
            pair_init[i] = triangles_src_[locally_matched_pairs_[i_grp].first][i].vertex;
        for (int i = 0; i < 3; ++i) // target triangle
            pair_init[i + 3] = triangles_tgt_[locally_matched_pairs_[i_grp].second][i].vertex;

        // grow the group
        for (int i_pair = 0; i_pair < num_pairs; ++i_pair) {
            if (i_grp == i_pair) continue;
            // candidate triangle pair
            std::vector<int> pair_cand(6);
            for (int i = 0; i < 3; ++i) // source triangle
                pair_cand[i] = triangles_src_[locally_matched_pairs_[i_pair].first][i].vertex;
            for (int i = 0; i < 3; ++i) // target triangle
                pair_cand[i + 3] = triangles_tgt_[locally_matched_pairs_[i_pair].second][i].vertex;

            if (satisfyGlobalConsistency(pair_init, pair_cand))
                groups[i_grp].push_back(i_pair);
        }
        group_size[i_grp] = groups[i_grp].size(); // 保留符合条件的最多组数
    }

    // output the maximum group
    // if there are two maximum groups, just use the first one
    // std::max_element(group_size.begin(), group_size.end()) 返回的是一个迭代器，指向 group_size 中的最大元素。
    // 为了得到这个元素的 索引（即它在 group_size 中的位置），你需要计算这个迭代器与 group_size.begin() 之间的差值。差值就是最大元素的索引。
    int id_max_group = std::max_element(group_size.begin(), group_size.end()) - group_size.begin();
    // resize：这是容器的一个方法，用于改变容器的大小。
    // 如果新大小大于当前大小，容器会添加新元素（其值通常被初始化为默认值，如 0 对于数值类型，或者 nullptr 对于指针类型，
    // 具体取决于容器的类型和是否有提供初始化值）。
    // 如果新大小小于当前大小，则超出新大小范围的元素会被删除。
    globally_matched_pairs_.resize(group_size[id_max_group]);
    for (int i = 0; i < group_size[id_max_group]; ++i)
        globally_matched_pairs_[i] = groups[id_max_group][i];

    // output position matches
    std::set<std::pair<int, int>> set_matches;
    for (const auto& id_pair: globally_matched_pairs_) {
        const auto& id_tri_src = locally_matched_pairs_[id_pair].first; // .first == [0]
        const auto& id_tri_tgt = locally_matched_pairs_[id_pair].second;
        for (int i = 0; i < 3; ++i)
            // insert:在容器的指定位置插入一个元素或一组元素。
            set_matches.insert(std::make_pair(triangles_src_[id_tri_src][i].vertex,
                                              triangles_tgt_[id_tri_tgt][i].vertex));
    }
    stem_matches_.resize(set_matches.size());
    int idx = 0;
    for (const auto& match: set_matches) {
        stem_matches_[idx].index_query = match.first;
        stem_matches_[idx].index_match = match.second;
        idx++;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Matching::randomGlobalMatching() {
//	const auto& num_pairs = locally_matched_pairs_.size();
    const auto num_pairs = 10000;
    std::vector<std::vector<int>> groups(num_pairs);
    std::vector<int> group_size(num_pairs);

    // generate random numbers from 0 to num_pairs as a vector
    std::vector<int> random_numbers(locally_matched_pairs_.size());
    std::iota(random_numbers.begin(), random_numbers.end(), 0);
    std::shuffle(random_numbers.begin(), random_numbers.end(), std::mt19937(std::random_device()()));
    std::vector<int> samples(random_numbers.begin(), random_numbers.begin() + num_pairs);

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int i_grp = 0; i_grp < num_pairs; ++i_grp) {
//		groups[i_grp].push_back(i_grp);
        const auto& pairid_sample = samples[i_grp];
        groups[i_grp].push_back(pairid_sample);
        // initial triangle pair
        std::vector<int> pair_init(6);
        for (int i = 0; i < 3; ++i) // source triangle
            pair_init[i] = triangles_src_[locally_matched_pairs_[pairid_sample].first][i].vertex;
        for (int i = 0; i < 3; ++i) // target triangle
            pair_init[i + 3] = triangles_tgt_[locally_matched_pairs_[pairid_sample].second][i].vertex;

        // grow the group
        for (int i_pair = 0; i_pair < num_pairs; ++i_pair) {
            if (pairid_sample == i_pair) continue;
            // candidate triangle pair
            std::vector<int> pair_cand(6);
            for (int i = 0; i < 3; ++i) // source triangle
                pair_cand[i] = triangles_src_[locally_matched_pairs_[i_pair].first][i].vertex;
            for (int i = 0; i < 3; ++i) // target triangle
                pair_cand[i + 3] = triangles_tgt_[locally_matched_pairs_[i_pair].second][i].vertex;

            if (satisfyGlobalConsistency(pair_init, pair_cand))
                groups[i_grp].push_back(i_pair);
        }
        group_size[i_grp] = groups[i_grp].size();
    }

    // output the maximum group
    // if there are two maximum groups, just use the first one
    int id_max_group = std::max_element(group_size.begin(), group_size.end()) - group_size.begin();
    globally_matched_pairs_.resize(group_size[id_max_group]);
    for (int i = 0; i < group_size[id_max_group]; ++i)
        globally_matched_pairs_[i] = groups[id_max_group][i];

    // output position matches
    std::set<std::pair<int, int>> set_matches;
    for (const auto& id_pair: globally_matched_pairs_) {
        const auto& id_tri_src = locally_matched_pairs_[id_pair].first;
        const auto& id_tri_tgt = locally_matched_pairs_[id_pair].second;
        for (int i = 0; i < 3; ++i)
            set_matches.insert(std::make_pair(triangles_src_[id_tri_src][i].vertex,
                                              triangles_tgt_[id_tri_tgt][i].vertex));
    }
    stem_matches_.resize(set_matches.size());
    int idx = 0;
    for (const auto& match: set_matches) {
        stem_matches_[idx].index_query = match.first;
        stem_matches_[idx].index_match = match.second;
        idx++;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Matching::estimateTransformation(Eigen::Matrix4f& transform) {
    constructTriangles(stem_positions_src_, triangles_src_);
    constructTriangles(stem_positions_tgt_, triangles_tgt_);

    localMatching();

    globalMatching();
    //randomGlobalMatching();

    // the 4-DoF solution
    // calculate 2D transformation matrix
    pcl::registration::TransformationEstimation2D<Point3D, Point3D, float> te;
    te.estimateRigidTransformation(
            *stem_positions_src_, *stem_positions_tgt_, stem_matches_, transform);

    // calculate height translation
    const auto& count_matches = stem_matches_.size();
    float z_src = 0.0, z_tgt = 0.0;
    for (int i = 0; i < count_matches; ++i) {
        z_src += stem_positions_src_->points[stem_matches_[i].index_query].z;
        z_tgt += stem_positions_tgt_->points[stem_matches_[i].index_match].z;
    }
    z_src = z_src / float(count_matches);
    z_tgt = z_tgt / float(count_matches);
    float trl_z = z_tgt - z_src;
    transform(2, 3) = trl_z;
}
