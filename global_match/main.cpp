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

#include <iomanip>
// PCL
#include <pcl/io/ply_io.h>
#include <pcl/console/parse.h>
// GlobalMatch
#include "stem_mapping.h"
#include "stem_matching.h"

using namespace pcl::console;

void
writeTransformationMatrix(const std::string& output_file, // const:保证了参数在函数内部不会被修改
                          const Eigen::Matrix4f& transformation_matrix) {
    std::ofstream outFile(output_file);
    if (outFile.is_open()) {
        outFile << std::setprecision(8) << std::fixed << transformation_matrix; // write matrix to file
        outFile.close();
    } else {
        // handle error or throw exception
        std::cerr << "Unable to open file: " + output_file << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char** argv) {                         // console parse, char** argv = char* argv[]
// Support multi-threaded compilation and linking
#if defined(_OPENMP)                                  // "D:\CxxProjects\GlobalMatch\GlobalMatch\code\CMakeLists.txt"
    print_info("[PARALLEL PROCESSING USING ");
    print_value("%d", omp_get_max_threads());
    print_info(" THREADS] \n\n");
#else
    print_info("[NON-PARALLEL PROCESSING] \n\n");
#endif
    // from consle prompt
    // argv[0] is the name of the .exe
    std::string filename_source = argv[1];
    std::string filename_target = argv[2];
    std::string filename_result = argv[3];

    double tic, toc, time_val = 0.0;

    //============================= Load Data =============================
    std::cout << "1. LOADING RAW POINT CLOUDS:" << std::endl;
    std::cout << "Source filename: " << filename_source << std::endl;
    std::cout << "Target filename: " << filename_target << std::endl;
    std::cout << "Result filename: " << filename_result << std::endl;

    Cloud3D::Ptr cloud_src(new Cloud3D), cloud_tgt(new Cloud3D); // "Cloud3D" defined in common.h; cloud_src: pointer; new Cloud3D: instance object
    tic = omp_get_wtime(); // _OPENMP
    pcl::io::loadPLYFile<Point3D>(filename_source, *cloud_src);
    toc = omp_get_wtime();
    print_info("  (1) ");
    print_value("%d", cloud_src->size()); // %d: int
    print_info(" source points in ");
    print_value("%f", toc - tic); // %f: float
    print_info(" s.\n");
    tic = omp_get_wtime();
    pcl::io::loadPLYFile<Point3D>(filename_target, *cloud_tgt);
    toc = omp_get_wtime();
    print_info("  (2) ");
    print_value("%d", cloud_tgt->size());
    print_info(" target points in ");
    print_value("%f", toc - tic);
    print_info(" s.\n");

    //======================= Pairwise Registration =======================
    std::cout << "2. MAPPING STEMS:" << std::endl;
    Cloud3D::Ptr cloud_pos_src(new Cloud3D), cloud_pos_tgt(new Cloud3D);
    tic = omp_get_wtime();
    Mapping mapping;
    // mapping.setInputCloud(std::shared_ptr<Cloud3D>(new Cloud3D(*cloud_src)));
    mapping.setInputCloud(cloud_src->makeShared()); // Will not change the value of cloud_strc or the object it points to
    mapping.extract(cloud_pos_src);
    toc = omp_get_wtime();
    time_val += (toc - tic);
    print_info("  (1) ");
    print_value("%d", cloud_pos_src->size());
    print_info(" source stem positions in ");
    print_value("%f", toc - tic);
    print_info(" s.\n");
    tic = omp_get_wtime();
    mapping.setInputCloud(cloud_tgt->makeShared());
    mapping.extract(cloud_pos_tgt);
    toc = omp_get_wtime();
    time_val += (toc - tic);
    print_info("  (2) ");
    print_value("%d", cloud_pos_tgt->size());
    print_info(" target stem positions in ");
    print_value("%f", toc - tic);
    print_info(" s.\n");

    std::cout << "3. MATCHING STEMS:" << std::endl;
    tic = omp_get_wtime();
    Eigen::Matrix4f mat_crs;
    Matching matching;
    matching.setPairwiseStemPositions(cloud_pos_src,
                                      cloud_pos_tgt);
    matching.estimateTransformation(mat_crs);
    toc = omp_get_wtime();
    time_val += (toc - tic);
    print_info("   ");
    print_value("%d", matching.getNumberOfMatches());
    print_info(" pairs of correspondences in ");
    print_value("%f", toc - tic);
    print_info(" s.\n");

    print_info("====> [Total running time] ");
    print_value("%f", time_val);
    print_info(" s.\n");

    //============================ Output Matrix ==========================
    writeTransformationMatrix(filename_result, mat_crs);

    return 0;
}

