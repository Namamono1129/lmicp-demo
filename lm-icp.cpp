#pragma once

#include "stdafx.h"
#include "projection_parameters.h"
#include <iostream>
#include "eigen/Eigen/Dense"

ProjectionParameters FitByLmIcp(ProjectionParameters initialParams, int N_d, Eigen::VectorXf pointsToFit, Eigen::VectorXf groundTruthPoints)
{
	const int N_p = 12;
	const float epsilon = 0.1;
	const int max_iterations = 100;
	float lambda = 0.001;
	float max_lambda = 10.0f;

	Eigen::VectorXf vec_a = initialParams.GetParamsAsVector();	// set a = a_0
	auto params = initialParams.createFromVector(vec_a);
	int k = 0;
	bool nan_detected = false;
	while (k < max_iterations && lambda < max_lambda && !nan_detected)
	{
		// Compute e_k (Error Vector)
		Eigen::VectorXf projectedPoints = params.Project(N_d, pointsToFit);

		Eigen::VectorXf X(N_d);
		Eigen::VectorXf Y(N_d);
		Eigen::VectorXf vec_e_k(N_d);
		for (int i = 0; i < N_d; i++)
		{
			X(i) = groundTruthPoints(i * 2 + 0) - projectedPoints(i * 2 + 0);
			Y(i) = groundTruthPoints(i * 2 + 1) - projectedPoints(i * 2 + 1);
			vec_e_k(i) = sqrt(X(i) * X(i) + Y(i) * Y(i));
		}

		// Compute Jacobian Matrix
		// For experiment and limiting parameters into o_x & o_y, now J is {Nd~2}

		Eigen::MatrixXf J(N_d, N_p);

		for (int i = 0; i < N_d; i++)
		{
			float sqrtPowXPlusPowY = sqrt(X(i)*X(i) + Y(i)*Y(i));
			float x_i = pointsToFit(i * 3 + 0);
			float y_i = pointsToFit(i * 3 + 1);
			float z_i = pointsToFit(i * 3 + 2);
			float r_11 = params.GetParamsAsVector()(0);
			float r_12 = params.GetParamsAsVector()(1);
			float r_13 = params.GetParamsAsVector()(2);
			float r_21 = params.GetParamsAsVector()(3);
			float r_22 = params.GetParamsAsVector()(4);
			float r_23 = params.GetParamsAsVector()(5);
			float f = params.GetParamsAsVector()(6);
			float t_x = params.GetParamsAsVector()(7);
			float t_y = params.GetParamsAsVector()(8);
			// r1-r6
			J(i, 0) = -X(i) / sqrtPowXPlusPowY * f / z_i * x_i;
			J(i, 1) = -X(i) / sqrtPowXPlusPowY * f / z_i * y_i;
			J(i, 2) = -X(i) * f / sqrtPowXPlusPowY;
			J(i, 3) = Y(i) / sqrtPowXPlusPowY * f / z_i * x_i;
			J(i, 4) = Y(i) / sqrtPowXPlusPowY * f / z_i * y_i;
			J(i, 5) = Y(i) * f / sqrtPowXPlusPowY;
			// f
			J(i, 6) = (-X(i)*((r_11*x_i + r_12 * y_i + r_13 * z_i + t_x) / z_i) + Y(i)*((r_21*x_i + r_22 * y_i + r_23 * z_i + t_y) / z_i)) / sqrtPowXPlusPowY;
			// tau
			J(i, 7) = -X(i) * f / sqrtPowXPlusPowY / z_i;	// ÝEi/ÝƒÑ_x
			J(i, 8) = Y(i) * f / sqrtPowXPlusPowY / z_i;		// ÝEi/ÝƒÑ_y
			J(i, 9) = 0;															// ÝEi/ÝƒÑ_z
																					// o
			J(i, 10) = -X(i) / sqrtPowXPlusPowY;	// ÝEi/Ýo_x
			J(i, 11) = -Y(i) / sqrtPowXPlusPowY;	// ÝEi/Ýo_y
		}
		Eigen::MatrixXf JT = J.transpose();

		while (lambda < max_lambda)
		{
			Eigen::MatrixXf lambdaI = Eigen::MatrixXf::Identity(N_p, N_p) * lambda;
			Eigen::VectorXf vec_a_k = vec_a - (JT * J + lambdaI).inverse() * JT * vec_e_k;
			// Compute New Errors
			Eigen::VectorXf vec_e_k_candidate(N_d);
			auto newParams = params.createFromVector(vec_a_k);
			Eigen::VectorXf newProjectedPoints = newParams.Project(N_d, pointsToFit);
			for (int i = 0; i < N_d; i++)
			{
				X(i) = groundTruthPoints(i * 2 + 0) - newProjectedPoints(i * 2 + 0);
				Y(i) = groundTruthPoints(i * 2 + 1) - newProjectedPoints(i * 2 + 1);
				vec_e_k_candidate(i) = sqrt(X(i) * X(i) + Y(i) * Y(i));
				if (isnan(vec_e_k_candidate(i)))
				{
					nan_detected = true;
				}
			}
			std::cout << "lambda=" << lambda << std::endl;
			std::cout << "e_k:" << std::endl;
			std::cout << vec_e_k << std::endl;
			std::cout << "new e_k:" << std::endl;
			std::cout << vec_e_k_candidate << std::endl;
			std::cout << "new a:" << std::endl;
			std::cout << vec_a_k << std::endl;

			if (nan_detected)
			{
				std::cout << "NaN detected: break" << std::endl;
				break;
			}

			if (vec_e_k.norm() > vec_e_k_candidate.norm())
			{
				std::cout << "break!" << std::endl << std::endl;
				params = params.createFromVector(vec_a_k);
				vec_a = params.GetParamsAsVector();
				break;
			}
			lambda = lambda * 2;
		}
		k++;
	}

	return params;
}
