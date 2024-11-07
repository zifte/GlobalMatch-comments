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

#include "stem_mapping.h"
 // pcl
#include <pcl/common/angles.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
// remove the leaf-size check
#include "voxel_grid_fix.h"
#include "octree_extract_clusters.h"
// csf
#include "CSF.h" // https://github.com/jianboqi/CSF
// cc
#include "jacobi.h"
// ICRA 2015 octree
#include "octree_unibn.hpp"
// vtk
#include <vtkOBBTree.h>

// transformCylinderAxisToTwoPoints():
// Converts a point-slope equation representing the direction of a cylinder's axis into two specific 3D points
// on that axis. These points are often used to define or represent the bottom and top locations of a cylinder.

// 当你写 double pt_bottom[3] 作为函数参数时，它并不会像你所期望的那样直接传递数组本身给函数。
// 相反，数组名 pt_bottom 在这种情况下会退化为指向数组第一个元素的指针（即 double*）。
// 这意味着函数接收的是一个指向 double 的指针，而不是数组的完整副本，但也不是数组的引用。
// 由于数组名退化为指针，函数内部无法知道原始数组的大小（除非这个信息以其他方式传递给函数，比如作为额外的参数）。
// 因此，使用 double pt_bottom[3] 作为函数参数时，你失去了数组大小的信息，
// 并且无法像使用数组引用那样在编译时检查数组大小是否匹配。
void
transformCylinderAxisToTwoPoints(const Coefficient& coeff,  // input
	const double z_bottom,     // input
	const double z_top,        // input
	double(&pt_bottom)[3],   // output
	double(&pt_top)[3]) {    // output
	// transform the point-slope line coefficients to two 3D points
	const auto& x0 = coeff.values[0]; // 表示x0是一个引用，而不是值的拷贝,以避免不必要的拷贝开销.
	const auto& y0 = coeff.values[1];
	const auto& z0 = coeff.values[2];
	const auto& dx = coeff.values[3];
	const auto& dy = coeff.values[4];
	const auto& dz = coeff.values[5];
	pt_bottom[0] = (z_bottom - z0) * dx / dz + x0; // x_bottom
	pt_bottom[1] = (z_bottom - z0) * dy / dz + y0; // y_bottom
	pt_bottom[2] = z_bottom;
	pt_top[0] = (z_top - z0) * dx / dz + x0; // x_top
	pt_top[1] = (z_top - z0) * dy / dz + y0; // y_top
	pt_top[2] = z_top;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
fitOneStem(const Cloud3D::Ptr& cloud_stem, // input
	int min_pts_per_cylinder, // input
	pcl::ModelCoefficients::Ptr& coefficients_cylinder) { // output
	// estimate normals: 法线提供了点云中每个点的表面方向信息。
	// 在拟合过程中，法线的方向可以帮助算法更准确地识别出符合圆柱体模型的点。
	// 法线可用于约束拟合过程，确保拟合的圆柱体与点云的表面形状一致。
	// 例如，通过设置法线与圆柱体轴向的最大夹角，确保拟合结果的合理性。
	CloudNormal::Ptr normals(new CloudNormal);
	pcl::search::KdTree<Point3D>::Ptr tree(new pcl::search::KdTree<Point3D>);
	tree->setInputCloud(cloud_stem);
	pcl::NormalEstimation<Point3D, pcl::Normal> ne;
	ne.setSearchMethod(tree);
	ne.setInputCloud(cloud_stem);
	ne.setKSearch(50);
	ne.compute(*normals);

	// cylinder fitting
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
	pcl::SACSegmentationFromNormals<Point3D, pcl::Normal> seg;
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_CYLINDER);
	seg.setAxis(Eigen::Vector3f(0.0, 0.0, 1.0));
	seg.setEpsAngle(pcl::deg2rad(50.0f));
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setNormalDistanceWeight(0.1);
	seg.setMaxIterations(10000);
	seg.setDistanceThreshold(0.03);
	seg.setRadiusLimits(0.05, 0.50);
	seg.setInputCloud(cloud_stem);
	seg.setInputNormals(normals);
	seg.segment(*inliers, *coefficients_cylinder);

	// check the quality of the cylinder
	if (inliers->indices.size() >= min_pts_per_cylinder)
		return true;
	else
		return false;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 下层植被通常指的是地面以上但低于树冠层的植被，比如灌木丛、低矮的树木等。
void
Mapping::extractUnderstory() {
	std::vector<int> indices_extracted, indices_remained;
	CSF csf;
	csf.params.bSloopSmooth = true; // "cloth_nodes.txt" will be created
	csf.params.dist_max = 3; // 将距离 DTM 0.2–3 米内的点标记为林下层。
	csf.params.dist_min = 0.2;
	csf.params.cloth_resolution = 0.5;
	csf.params.interations = 500;
	csf.params.rigidness = 2;
	csf.setPointCloud(cloud_input_);
	csf.do_filtering(indices_extracted, indices_remained, true, mesh_ground_);
	indices_understory_->indices = indices_extracted;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Mapping::extractStemPoints() {
	// downsampling
	Cloud3D::Ptr cloud_understory(new Cloud3D);
	pcl::VoxelGrid<Point3D> vg;
	vg.setInputCloud(cloud_input_);
	vg.setIndices(indices_understory_);
	vg.setLeafSize(leaf_size_, leaf_size_, leaf_size_);
	vg.filter(*cloud_understory);

	// remove understory shrub based on pointwise verticality
	size_t cloud_size = cloud_understory->size(); // understory point cloud
	unibn::Octree<Point3D> octree;
	octree.initialize(*cloud_understory);
	// 这是一个二维向量（嵌套向量）。外层向量的每个元素都是一个内层向量，内层向量存储 uint32_t 类型的数据。
	// uint32_t 是无符号的 32 位整数，通常用于存储索引或标识符。
	// 通过构造函数初始化外层向量的大小为 cloud_size。
	// 每个内层向量的大小在此时是未定义的，但可以根据需要动态添加元素。
	// 如果想定义内层向量大小为100：std::vector<std::vector<uint32_t>> neighbors(cloud_size, std::vector<uint32_t>(100));
	std::vector<std::vector<uint32_t>> neighbors(cloud_size);

	// 这是一个一维向量，存储 float 类型的数据，用于表示浮点数值。
	// 通过构造函数初始化 verticality 的大小为 cloud_size，并将所有元素的初始值设置为 0.0
	std::vector<float> verticality(cloud_size, 0.0);
#if defined(_OPENMP)
#pragma omp parallel
#endif
	{ // #pragma omp parallel指令开启了一个并行区域。大括号内的代码将在多个线程中并行执行
		// query neighbors
#if defined(_OPENMP)
#pragma omp for
#endif
		for (int i = 0; i < cloud_size; ++i)
			octree.radiusNeighbors<unibn::L2Distance<Point3D>>(
				cloud_understory->points[i], search_radius_, neighbors[i]);

		// calculate point-wise verticality
#if defined(_OPENMP)
#pragma omp for
#endif

		for (int i = 0; i < neighbors.size(); ++i) {
			// need at least 3 points
			// (noted that i-th point's neighbors include itself)
			if (neighbors[i].size() < 4)
				continue;

			// calculate the gravity center of i-th point's neighbors
			const auto& count = neighbors[i].size();
			// {} 是正确的方式来初始化数组。
			// [] 用于声明数组或访问数组元素，但不能用于初始化。
			float gravity_center[3] = { 0, 0, 0 };
			for (auto& neighbor : neighbors[i]) {
				gravity_center[0] += cloud_understory->points[neighbor].x;
				gravity_center[1] += cloud_understory->points[neighbor].y;
				gravity_center[2] += cloud_understory->points[neighbor].z;
			}
			gravity_center[0] /= count;
			gravity_center[1] /= count;
			gravity_center[2] /= count;

			// build a covariance matrix
			float mxx = 0.0, myy = 0.0, mzz = 0.0, mxy = 0.0, mxz = 0.0, myz = 0.0;
			for (auto& neighbor : neighbors[i]) {
				float dx = cloud_understory->points[neighbor].x - gravity_center[0];
				float dy = cloud_understory->points[neighbor].y - gravity_center[1];
				float dz = cloud_understory->points[neighbor].z - gravity_center[2];
				mxx += dx * dx;
				myy += dy * dy;
				mzz += dz * dz;
				mxy += dx * dy;
				mxz += dx * dz;
				myz += dy * dz;
			}
			CCCoreLib::SquareMatrixf mat_cov(3);
			mat_cov.m_values[0][0] = mxx / count;
			mat_cov.m_values[1][1] = myy / count;
			mat_cov.m_values[2][2] = mzz / count;
			mat_cov.m_values[1][0] = mat_cov.m_values[0][1] = mxy / count;
			mat_cov.m_values[2][0] = mat_cov.m_values[0][2] = mxz / count;
			mat_cov.m_values[2][1] = mat_cov.m_values[1][2] = myz / count;

			CCCoreLib::SquareMatrixf eigen_vectors;
			std::vector<float> eigen_values;
			if (!CCCoreLib::Jacobi<float>::ComputeEigenValuesAndVectors(
				mat_cov, eigen_vectors, eigen_values, true))
				// failed to compute the eigen values
				continue;

			// sort the eigenvectors in decreasing order of their associated eigenvalues
			// 这意味着最大的特征值对应的特征向量将排在前面，这通常代表着数据在该方向上的变化性或重要性。
			// 特征值: 反映了在特定方向上数据的分散程度。特征值越大，表示在该方向上数据的变异性越大。
		    // 特征向量: 表示数据变化的方向。排序后，前面的特征向量对应的特征值更大，意味着在这些方向上，数据的变化更显著。
			CCCoreLib::Jacobi<float>::SortEigenValuesAndVectors(
				eigen_vectors, eigen_values);
			// 定义了一个三维向量 z，表示垂直方向（Z轴方向）。
			CCVector3f z(0, 0, 1);
			// 初始化一个新的向量 e3，并将其设置为 z 的副本。
			CCVector3f e3(z);
			// 从已排序的特征向量 eigen_vectors 中提取第三个特征向量，并将其存储在 e3.u 中。
			CCCoreLib::Jacobi<float>::GetEigenVector(eigen_vectors, 2, e3.u);
			// 计算垂直向量 z（即 Z 轴方向）和特征向量 e3 的点积。
			// 点积用于衡量两个向量之间的角度关系。
			// 点积的绝对值，结果在 0 到 1 之间,
			// 值接近 1 表示两个向量接近同一方向，值接近 0 表示接近垂直。
			verticality[i] = 1.0f - std::abs(z.dot(e3));
		}
	}
	for (int i = 0; i < verticality.size(); ++i) {
		if (verticality[i] >= verticality_threshold_)
			cloud_stems_->push_back(cloud_understory->points[i]);
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Mapping::extractTreePositions(Cloud3D::Ptr& tree_positions) {
	// Euclidean cluster extraction
	std::vector<pcl::PointIndices> ptid_clusters;
	OctreeEuclideanClusterExtraction<Point3D> oec;
	oec.setClusterTolerance(min_dist_between_stems_);
	oec.setMinClusterSize(min_pts_per_cluster_);
	oec.setInputCloud(cloud_stems_);
	oec.extract(ptid_clusters);

	// Cylinder fitting
	CoefficientVector coefficients(ptid_clusters.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
	for (int i = 0; i < ptid_clusters.size(); ++i) {
		Cloud3D::Ptr cloud_cluster(new Cloud3D);
		for (int& pt_index : ptid_clusters[i].indices)
			cloud_cluster->push_back((*cloud_stems_)[pt_index]); // point coordinates of a single cluster

		pcl::ModelCoefficients::Ptr coefficients_cylinder(new pcl::ModelCoefficients);
		if (fitOneStem(cloud_cluster, min_pts_per_cluster_, coefficients_cylinder))
			coefficients[i] = *coefficients_cylinder;
	}

	// Extract candidate positions
	// calculate two given z values
	pcl::PointCloud<pcl::PointXYZL>::Ptr candidate_positions(new pcl::PointCloud<pcl::PointXYZL>);
	double bbox[6]; // x_min, x_max, y_min, y_max, z_min, z_max
	mesh_ground_->GetBounds(bbox);
	double z_bottom = bbox[4] - 50.0; // z_min - rough guess， 目的：两顶点连线要穿过DEM，便于计算圆柱轴线与DEM的交点
	double z_top = bbox[5] + 50.0; // z_max + rough guess
	// vtkOBBTree 是 VTK 中的一个类，表示一个有序包围盒树（Oriented Bounding Box Tree）。
	// 可以用vtkOBBTree和直线、三角形甚至是另一个vtkOBBTree做相交检测、运算，碰撞检测
	// 可以用来，获得直线与多边形数据的交点,用来求最近点,等
	vtkSmartPointer<vtkOBBTree> obbtree = vtkSmartPointer<vtkOBBTree>::New();
	obbtree->SetDataSet(mesh_ground_);
	obbtree->BuildLocator();
	for (int id_cluster = 0; id_cluster < ptid_clusters.size(); ++id_cluster) {
		if (coefficients[id_cluster].values.empty())
			continue;
		// transform the point-slope line coefficients to two 3D points
		double pt_bottom[3], pt_top[3];
		transformCylinderAxisToTwoPoints(coefficients[id_cluster], z_bottom, z_top, pt_bottom, pt_top);

		// 初始化一个用于存储空间中多点的容器。
		// 在具体应用中，intersections 可能会用于记录如相交点、模型顶点、路径节点等点集信息。
		vtkSmartPointer<vtkPoints> intersections = vtkSmartPointer<vtkPoints>::New();
		// 计算圆柱轴线与DEM的交点：obbtree中已经存储了DEM模型信息
		// pt_bottom 和 pt_top：这两个参数是定义线段的起点和终点的 3D 坐标。
		// intersections：这是一个 vtkPoints 类型的对象，用于存储线段与模型的所有交点。
		// nullptr：指定即交点所属的模型单元 ID。这里使用 nullptr 表示忽略该信息。
		obbtree->IntersectWithLine(pt_bottom, pt_top, intersections, nullptr);

		if (intersections->GetNumberOfPoints() == 1) { // 只有一个交点
			double intersection[3];
			intersections->GetPoint(0, intersection);
			pcl::PointXYZL position_with_cluster_id; // L: label
			position_with_cluster_id.x = static_cast<float>(intersection[0]);
			position_with_cluster_id.y = static_cast<float>(intersection[1]);
			position_with_cluster_id.z = static_cast<float>(intersection[2]);
			position_with_cluster_id.label = id_cluster;
			candidate_positions->push_back(position_with_cluster_id);
		}
	}

	// Optimize stem positions (to avoid generating several positions representing the same stem)
	// 这是一个标准库中的向量（std::vector），其中每个元素的类型是 pcl::PointIndices。
	std::vector<pcl::PointIndices> stemid_clusters;
	OctreeEuclideanClusterExtraction<pcl::PointXYZL> oec_stem;
	oec_stem.setClusterTolerance(min_dist_between_stems_);
	oec_stem.setInputCloud(candidate_positions);
	oec_stem.extract(stemid_clusters);

	for (auto& stemid_cluster : stemid_clusters) {
		if (stemid_cluster.indices.size() == 1) { // qualified stem positions
			const auto& stemid = stemid_cluster.indices[0];
			Point3D position;
			position.x = candidate_positions->points[stemid].x;
			position.y = candidate_positions->points[stemid].y;
			position.z = candidate_positions->points[stemid].z;
			tree_positions->push_back(position);
		}
		else { // unqualified stem positions (that are too close to each other)
			Cloud3D::Ptr merged_cluster(new Cloud3D);
			for (auto& stemid : stemid_cluster.indices) {
				const auto& cluster_id = candidate_positions->points[stemid].label;
				for (auto& ptid : ptid_clusters[cluster_id].indices)
					merged_cluster->push_back((*cloud_stems_)[ptid]);
			}
			// re-estimate a stem position from the merged cluster
			pcl::ModelCoefficients::Ptr coefficients_cylinder(new pcl::ModelCoefficients);
			if (fitOneStem(merged_cluster, min_pts_per_cluster_, coefficients_cylinder)) {
				// transform the point-slope line coefficients to two 3D points
				double pt_bottom[3], pt_top[3];
				transformCylinderAxisToTwoPoints(
					*coefficients_cylinder, z_bottom, z_top, pt_bottom, pt_top);
				vtkSmartPointer<vtkPoints> intersections = vtkSmartPointer<vtkPoints>::New();
				obbtree->IntersectWithLine(pt_bottom, pt_top, intersections, nullptr);
				if (intersections->GetNumberOfPoints() == 1) {
					double intersection[3];
					intersections->GetPoint(0, intersection);
					Point3D position;
					position.x = static_cast<float>(intersection[0]);
					position.y = static_cast<float>(intersection[1]);
					position.z = static_cast<float>(intersection[2]);
					tree_positions->push_back(position);
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Mapping::extract(Cloud3D::Ptr& tree_positions) {
	if (!cloud_input_) {
		PCL_WARN("[Mapping::extract] No input dataset given!\n");
		tree_positions->width = tree_positions->height = 0;
		tree_positions->points.clear();
		return;
	}

	// initialize
	tree_positions->clear();
	indices_understory_->indices.clear();
	cloud_stems_->clear();
	// mesh_ground_ = vtkSmartPointer<vtkPolyData>::New();
	mesh_ground_->Initialize();

	Indices::Ptr indices_understory(new Indices);
	extractUnderstory();
	extractStemPoints();
	extractTreePositions(tree_positions);
}
