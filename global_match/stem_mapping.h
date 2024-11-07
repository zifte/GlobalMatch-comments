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

#ifndef GLOBALMATCH_STEM_MAPPING_H
#define GLOBALMATCH_STEM_MAPPING_H

#include "common.h"
#include <pcl/features/normal_3d.h>
#include <vtkSmartPointer.h> // Visualization Toolkit
#include <vtkPolyData.h>

typedef pcl::PointIndices Indices; // Stores the index of the points
typedef pcl::ModelCoefficients Coefficient; // Stores model coefficients such as plane, sphere, cylinder
typedef std::vector<Coefficient> CoefficientVector; // Stores multiple pcl::ModelCoefficients objects

class Mapping { // class
public:
    // Mapping(): Constructor; used to set default parameter values
    Mapping() : leaf_size_(0.01), // defined in private
                search_radius_(0.10), // defined in private
                verticality_threshold_(0.90), // defined in private
                min_pts_per_cluster_(200), // defined in private
                min_dist_between_stems_(2 * search_radius_), // defined in private
                mesh_ground_(vtkSmartPointer<vtkPolyData>::New()), // defined in private, Automatic memory management
                indices_understory_(new Indices), // defined in private,Manual memory management
                cloud_stems_(new Cloud3D) { // defined in private,Manual memory management
    }

    // ~Mapping(): Destructor: A special member function that performs cleanup work at the end of an object's lifetime, 
    //             such as releasing resources that the object allocated during its lifetime. 
    //             The name of the Destructor is the same as the class name, but with a tilde (~) in front of it.
    // virtual:    Allow derived classes to override functions of the same name in parent classes.
    // default:    Generates a default destructor implementation for Mapping().

    virtual ~Mapping() = default;

    // inline: Embed function code directly into the place where the function is called

    // const Cloud3D pointCloud; 
    // Cloud3D::ConstPtr cloudInput = std::make_shared<const Cloud3D>(pointCloud); 
    // 这里const表示pointCloud是不可以被修改的，
    // ConstPtr 表示不能通过cloudInput 来修改pointCloud
    // const: 确保Cloud3D对象本身是只读的。
    // ConstPtr: 作为智能指针类型，保证了无法通过这个指针对其指向的对象进行修改。
    inline void 
    setInputCloud(const Cloud3D::ConstPtr& cloud_input) {
        cloud_input_ = cloud_input;
    }

    inline void //Member function of Mapping(); use inline keyword to reduce function call overhead
    setLeafSize(float leaf_size) {
        leaf_size_ = leaf_size;
    }

    inline void //Member function of Mapping(); use inline keyword to reduce function call overhead
    setSearchRadius(float search_radius) {
        search_radius_ = search_radius;
    }

    inline void //Member function of Mapping(); use inline keyword to reduce function call overhead
    setVerticalityThreshold(float verticality_threshold) {
        verticality_threshold_ = verticality_threshold;
    }

    inline void //Member function of Mapping(); use inline keyword to reduce function call overhead
    setMinClusterSize(int min_cluster_size) {
        min_pts_per_cluster_ = min_cluster_size;
    }

    // TODO: add more setters and getters

    void extract(Cloud3D::Ptr& tree_positions);  //Member function of Mapping(), defined in stem_mapping.cpp

private:

    void extractUnderstory(); // defined in stem_mapping.cpp

    void extractStemPoints(); // defined in stem_mapping.cpp

    void extractTreePositions(Cloud3D::Ptr& tree_positions); // defined in stem_mapping.cpp

    Cloud3D::ConstPtr cloud_input_;
    Indices::Ptr indices_understory_;
    Cloud3D::Ptr cloud_stems_;
    vtkSmartPointer<vtkPolyData> mesh_ground_;

    // Subsampling
    float leaf_size_;
    // Verticality-based filtering
    float search_radius_;
    float verticality_threshold_;
    // Euclidean clustering
    int min_pts_per_cluster_;
    float min_dist_between_stems_;
};


#endif //GLOBALMATCH_STEM_MAPPING_H
