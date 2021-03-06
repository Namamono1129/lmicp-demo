// LMICP.cpp : アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"
#include "projection_parameters.h"
#include "lm-icp.h"
#include "sample/sample_pointclouds.h"
#include "sample/sample_parameters.h"
#include <iostream>
#include "eigen/Eigen/Dense"
#include "opencv/cv.hpp"

int main()
{
	const int N_d = 4;		// Number of data
	const int N_p = 12;		// Number of camera parameters

	/*** Load sample data points ***/
	Eigen::VectorXf pointsToProject = GetSamplePointsToProject(N_d);
	Eigen::VectorXf vec_xI = GetSampleGroundTruthPointCloud(N_d);

	/*** Load sample(default) projection parameters ***/
	ProjectionParameters originalParams = GetSampleProjectionParameters();
	ProjectionParameters paramsToFit = GetSampleProjectionParameters();

	/*** Fitting ***/
	ProjectionParameters fittedParams = FitByLmIcp(paramsToFit, N_d, pointsToProject, vec_xI);

	/*** Draw Points ***/
	cv::Mat im = cv::Mat::zeros(480, 640, CV_8UC3);

	Eigen::VectorXf originalProjectedPoints = originalParams.Project(N_d, pointsToProject);
	Eigen::VectorXf icpProjectedPoints = fittedParams.Project(N_d, pointsToProject);

	for (int i = 0; i < N_d; i++)
	{
		// Original Points with default projection parameters(Green)
		cv::circle(im, cv::Point(originalProjectedPoints[i * 2 + 0], originalProjectedPoints[i * 2 + 1]), 4, cv::Scalar(0, 255, 0), 3, CV_AA);
		// Ground-truth 2D points(blue)
		cv::circle(im, cv::Point(vec_xI(i * 2 + 0), vec_xI(i * 2 + 1)), 8, cv::Scalar(255, 0, 0), 3, CV_AA);
		// Projected 2D points fitted by LM-ICP
		cv::circle(im, cv::Point(icpProjectedPoints(i * 2 + 0), icpProjectedPoints(i * 2 + 1)), 4, cv::Scalar(0, 0, 255), 3, CV_AA);
	}
	cv::imshow("lm-icp demo", im);
	cv::waitKey(0);

    return 0;
}

