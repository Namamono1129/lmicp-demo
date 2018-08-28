#pragma once
#include "Eigen/Dense"

/*! Returns sample 2D ground-truth point cloud.
	@param[in]	N_d		Number of data points fitted in LM-ICP.
	@return				Vector of ground-truth 2D data points. {x_1, y_1, x_2, ...}
*/
Eigen::VectorXf GetSampleGroundTruthPointCloud(int N_d)
{
	float xI_1 = 480;
	float yI_1 = 272;
	float xI_2 = 576;
	float yI_2 = 332;
	float xI_3 = 320;
	float yI_3 = 332;
	float xI_4 = 220;
	float yI_4 = 280;
	Eigen::VectorXf vec_xI(N_d * 2);
	vec_xI(0) = xI_1;
	vec_xI(1) = yI_1;
	vec_xI(2) = xI_2;
	vec_xI(3) = yI_2;
	vec_xI(4) = xI_3;
	vec_xI(5) = yI_3;
	vec_xI(6) = xI_4;
	vec_xI(7) = yI_4;

	return vec_xI;
}

/*! Returns sample 3D point cloud to be projected.
	@param[in]	N_d		Number of data points fitted in LM-ICP.
	@return				Vector of 3D data points to be projected. {x_1, y_1, z_1, x_2, y_2, ...}
*/
Eigen::VectorXf GetSamplePointsToProject(int N_d)
{
	float x_1 = 0;
	float y_1 = -50;
	float z_1 = -50;
	float x_2 = -50;
	float y_2 = 50;
	float z_2 = -50;
	float x_3 = 50;
	float y_3 = 50;
	float z_3 = -50;
	float x_4 = 75;
	float y_4 = -50;
	float z_4 = -50;

	Eigen::VectorXf PointsToProject3D(N_d * 3);
	PointsToProject3D(0) = x_1;
	PointsToProject3D(1) = y_1;
	PointsToProject3D(2) = z_1;
	PointsToProject3D(3) = x_2;
	PointsToProject3D(4) = y_2;
	PointsToProject3D(5) = z_2;
	PointsToProject3D(6) = x_3;
	PointsToProject3D(7) = y_3;
	PointsToProject3D(8) = z_3;
	PointsToProject3D(9) = x_4;
	PointsToProject3D(10) = y_4;
	PointsToProject3D(11) = z_4;

	return PointsToProject3D;
}



