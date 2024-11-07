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

// ����д double pt_bottom[3] ��Ϊ��������ʱ��������������������������ֱ�Ӵ������鱾���������
// �෴�������� pt_bottom ����������»��˻�Ϊָ�������һ��Ԫ�ص�ָ�루�� double*����
// ����ζ�ź������յ���һ��ָ�� double ��ָ�룬�����������������������Ҳ������������á�
// �����������˻�Ϊָ�룬�����ڲ��޷�֪��ԭʼ����Ĵ�С�����������Ϣ��������ʽ���ݸ�������������Ϊ����Ĳ�������
// ��ˣ�ʹ�� double pt_bottom[3] ��Ϊ��������ʱ����ʧȥ�������С����Ϣ��
// �����޷���ʹ���������������ڱ���ʱ��������С�Ƿ�ƥ�䡣
void
transformCylinderAxisToTwoPoints(const Coefficient& coeff,  // input
	const double z_bottom,     // input
	const double z_top,        // input
	double(&pt_bottom)[3],   // output
	double(&pt_top)[3]) {    // output
	// transform the point-slope line coefficients to two 3D points
	const auto& x0 = coeff.values[0]; // ��ʾx0��һ�����ã�������ֵ�Ŀ���,�Ա��ⲻ��Ҫ�Ŀ�������.
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
	// estimate normals: �����ṩ�˵�����ÿ����ı��淽����Ϣ��
	// ����Ϲ����У����ߵķ�����԰����㷨��׼ȷ��ʶ�������Բ����ģ�͵ĵ㡣
	// ���߿�����Լ����Ϲ��̣�ȷ����ϵ�Բ��������Ƶı�����״һ�¡�
	// ���磬ͨ�����÷�����Բ������������нǣ�ȷ����Ͻ���ĺ����ԡ�
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
// �²�ֲ��ͨ��ָ���ǵ������ϵ��������ڲ��ֲ���������ľ�ԡ��Ͱ�����ľ�ȡ�
void
Mapping::extractUnderstory() {
	std::vector<int> indices_extracted, indices_remained;
	CSF csf;
	csf.params.bSloopSmooth = true; // "cloth_nodes.txt" will be created
	csf.params.dist_max = 3; // ������ DTM 0.2�C3 ���ڵĵ���Ϊ���²㡣
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
	// ����һ����ά������Ƕ�������������������ÿ��Ԫ�ض���һ���ڲ��������ڲ������洢 uint32_t ���͵����ݡ�
	// uint32_t ���޷��ŵ� 32 λ������ͨ�����ڴ洢�������ʶ����
	// ͨ�����캯����ʼ����������Ĵ�СΪ cloud_size��
	// ÿ���ڲ������Ĵ�С�ڴ�ʱ��δ����ģ������Ը�����Ҫ��̬���Ԫ�ء�
	// ����붨���ڲ�������СΪ100��std::vector<std::vector<uint32_t>> neighbors(cloud_size, std::vector<uint32_t>(100));
	std::vector<std::vector<uint32_t>> neighbors(cloud_size);

	// ����һ��һά�������洢 float ���͵����ݣ����ڱ�ʾ������ֵ��
	// ͨ�����캯����ʼ�� verticality �Ĵ�СΪ cloud_size����������Ԫ�صĳ�ʼֵ����Ϊ 0.0
	std::vector<float> verticality(cloud_size, 0.0);
#if defined(_OPENMP)
#pragma omp parallel
#endif
	{ // #pragma omp parallelָ�����һ���������򡣴������ڵĴ��뽫�ڶ���߳��в���ִ��
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
			// {} ����ȷ�ķ�ʽ����ʼ�����顣
			// [] ��������������������Ԫ�أ����������ڳ�ʼ����
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
			// ����ζ����������ֵ��Ӧ����������������ǰ�棬��ͨ�������������ڸ÷����ϵı仯�Ի���Ҫ�ԡ�
			// ����ֵ: ��ӳ�����ض����������ݵķ�ɢ�̶ȡ�����ֵԽ�󣬱�ʾ�ڸ÷��������ݵı�����Խ��
		    // ��������: ��ʾ���ݱ仯�ķ��������ǰ�������������Ӧ������ֵ������ζ������Щ�����ϣ����ݵı仯��������
			CCCoreLib::Jacobi<float>::SortEigenValuesAndVectors(
				eigen_vectors, eigen_values);
			// ������һ����ά���� z����ʾ��ֱ����Z�᷽�򣩡�
			CCVector3f z(0, 0, 1);
			// ��ʼ��һ���µ����� e3������������Ϊ z �ĸ�����
			CCVector3f e3(z);
			// ����������������� eigen_vectors ����ȡ����������������������洢�� e3.u �С�
			CCCoreLib::Jacobi<float>::GetEigenVector(eigen_vectors, 2, e3.u);
			// ���㴹ֱ���� z���� Z �᷽�򣩺��������� e3 �ĵ����
			// ������ں�����������֮��ĽǶȹ�ϵ��
			// ����ľ���ֵ������� 0 �� 1 ֮��,
			// ֵ�ӽ� 1 ��ʾ���������ӽ�ͬһ����ֵ�ӽ� 0 ��ʾ�ӽ���ֱ��
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
	double z_bottom = bbox[4] - 50.0; // z_min - rough guess�� Ŀ�ģ�����������Ҫ����DEM�����ڼ���Բ��������DEM�Ľ���
	double z_top = bbox[5] + 50.0; // z_max + rough guess
	// vtkOBBTree �� VTK �е�һ���࣬��ʾһ�������Χ������Oriented Bounding Box Tree����
	// ������vtkOBBTree��ֱ�ߡ���������������һ��vtkOBBTree���ཻ��⡢���㣬��ײ���
	// �������������ֱ�����������ݵĽ���,�����������,��
	vtkSmartPointer<vtkOBBTree> obbtree = vtkSmartPointer<vtkOBBTree>::New();
	obbtree->SetDataSet(mesh_ground_);
	obbtree->BuildLocator();
	for (int id_cluster = 0; id_cluster < ptid_clusters.size(); ++id_cluster) {
		if (coefficients[id_cluster].values.empty())
			continue;
		// transform the point-slope line coefficients to two 3D points
		double pt_bottom[3], pt_top[3];
		transformCylinderAxisToTwoPoints(coefficients[id_cluster], z_bottom, z_top, pt_bottom, pt_top);

		// ��ʼ��һ�����ڴ洢�ռ��ж���������
		// �ھ���Ӧ���У�intersections ���ܻ����ڼ�¼���ཻ�㡢ģ�Ͷ��㡢·���ڵ�ȵ㼯��Ϣ��
		vtkSmartPointer<vtkPoints> intersections = vtkSmartPointer<vtkPoints>::New();
		// ����Բ��������DEM�Ľ��㣺obbtree���Ѿ��洢��DEMģ����Ϣ
		// pt_bottom �� pt_top�������������Ƕ����߶ε������յ�� 3D ���ꡣ
		// intersections������һ�� vtkPoints ���͵Ķ������ڴ洢�߶���ģ�͵����н��㡣
		// nullptr��ָ��������������ģ�͵�Ԫ ID������ʹ�� nullptr ��ʾ���Ը���Ϣ��
		obbtree->IntersectWithLine(pt_bottom, pt_top, intersections, nullptr);

		if (intersections->GetNumberOfPoints() == 1) { // ֻ��һ������
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
	// ����һ����׼���е�������std::vector��������ÿ��Ԫ�ص������� pcl::PointIndices��
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
