#pragma once

#include "stdafx.h"
#include "projection_parameters.h"
#include "eigen/Eigen/Dense"

ProjectionParameters FitByLmIcp(ProjectionParameters initialParams, int N_d, Eigen::VectorXf pointsToFit, Eigen::VectorXf groundTruthPoints);