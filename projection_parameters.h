#pragma once
#include <iostream>
#include <cassert>
#include <Eigen/Dense>

class RotationMatrix
{
private:

public:
	float r_11;
	float r_12;
	float r_13;
	float r_21;
	float r_22;
	float r_23;

	RotationMatrix() {}

	RotationMatrix(float r_11, float r_12, float r_13, float r_21, float r_22, float r_23)
	{
		this->r_11 = r_11;
		this->r_12 = r_12;
		this->r_13 = r_13;
		this->r_21 = r_21;
		this->r_22 = r_22;
		this->r_23 = r_23;
	}

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> asBlockMatrix(int numData)
	{
		Eigen::MatrixXf outMatrix = Eigen::MatrixXf::Zero(numData * 3, numData * 3);
		for (int i = 0; i < numData; i++)
		{
			outMatrix(i * 3 + 0, i * 3 + 0) = this->r_11;
			outMatrix(i * 3 + 0, i * 3 + 1) = this->r_12;
			outMatrix(i * 3 + 0, i * 3 + 2) = this->r_13;
			outMatrix(i * 3 + 1, i * 3 + 0) = this->r_21;
			outMatrix(i * 3 + 1, i * 3 + 1) = this->r_22;
			outMatrix(i * 3 + 1, i * 3 + 2) = this->r_23;
			outMatrix(i * 3 + 2, i * 3 + 0) = 0;
			outMatrix(i * 3 + 2, i * 3 + 1) = 0;
			outMatrix(i * 3 + 2, i * 3 + 2) = 0;
		}
		return outMatrix;
	}
};

class CameraDisplacementVector
{
private:
	
public:
	float tau_x;
	float tau_y;
	float tau_z;

	CameraDisplacementVector(){}
	CameraDisplacementVector(float tau_x, float tau_y, float tau_z)
	{
		this->tau_x = tau_x;
		this->tau_y = tau_y;
		this->tau_z = tau_z;
	}

	Eigen::VectorXf asStackedVector(int numData)
	{
		Eigen::VectorXf outVector(3 * numData);
		for (int i = 0; i < numData; i++)
		{
			outVector(i * 3 + 0) = this->tau_x;
			outVector(i * 3 + 1) = this->tau_y;
			outVector(i * 3 + 2) = this->tau_z;
		}

		return outVector;
	}
};

class OpticalAxisVector
{
private:
	

public:
	float o_x;
	float o_y;

	OpticalAxisVector() {}
	OpticalAxisVector(float o_x, float o_y)
	{
		this->o_x = o_x;
		this->o_y = o_y;
	}

	Eigen::VectorXf asStackedVector(int numData)
	{
		Eigen::VectorXf outVector(2 * numData);
		for (int i = 0; i < numData; i++)
		{
			outVector(i * 2 + 0) = this->o_x;
			outVector(i * 2 + 1) = this->o_y;
		}
		return outVector;
	}
};

class ProjectionMatrix
{
private:
	
public:
	float f;

	ProjectionMatrix() {}

	ProjectionMatrix(float f)
	{
		this->f = f;
	}

	Eigen::MatrixXf asBlockMatrix(int numData, Eigen::VectorXf z)
	{
		Eigen::MatrixXf outMatrix = Eigen::MatrixXf::Zero(2 * numData, 3 * numData);
		for (int i = 0; i < numData; i++)
		{
			outMatrix(i * 2 + 0, i * 3 + 0) = this->f / z(i);
			outMatrix(i * 2 + 0, i * 3 + 1) = 0;
			outMatrix(i * 2 + 0, i * 3 + 2) = 0;
			outMatrix(i * 2 + 1, i * 3 + 0) = 0;
			outMatrix(i * 2 + 1, i * 3 + 1) = -(this->f / z(i));
			outMatrix(i * 2 + 1, i * 3 + 2) = 0;
		}
		return outMatrix;
	}
};

class ProjectionParameters
{
private:
	RotationMatrix rotationMatrix;
	float f;
	CameraDisplacementVector cameraDisplacementVector;
	OpticalAxisVector opticalAxisVector;
	ProjectionMatrix projectionMatrix;

public:
	ProjectionParameters(){}

	ProjectionParameters(
		RotationMatrix rotationMatrix,
		float f,
		CameraDisplacementVector cameraDisplacementVector,
		OpticalAxisVector opticalAxisVector)
	{
		this->rotationMatrix = rotationMatrix;
		this->f = f;
		this->cameraDisplacementVector = cameraDisplacementVector;
		this->opticalAxisVector = opticalAxisVector;
		this->projectionMatrix = ProjectionMatrix(f);
	}

	Eigen::VectorXf Project(int numData, Eigen::VectorXf pointsToProject)
	{
		assert(pointsToProject.size() == numData * 3);

		Eigen::VectorXf zs(numData);
		for (int i = 0; i < numData; i++)
		{
			zs(i) = pointsToProject(i * 3 + 2);
		}

		Eigen::VectorXf projectedPoints(numData * 2);
		Eigen::MatrixXf P = this->projectionMatrix.asBlockMatrix(numData, zs);
		Eigen::MatrixXf R = this->rotationMatrix.asBlockMatrix(numData);
		Eigen::VectorXf tau = this->cameraDisplacementVector.asStackedVector(numData);
		Eigen::VectorXf o = this->opticalAxisVector.asStackedVector(numData);
		projectedPoints = P * (R * pointsToProject + tau) + o;
		return projectedPoints;
	}

	void UpdateOpticalAxis(Eigen::Matrix2f new_vec_o)
	{
		this->opticalAxisVector = OpticalAxisVector(new_vec_o(0), new_vec_o(1));
	}


	Eigen::VectorXf GetParamsAsVector()
	{
		Eigen::VectorXf paramsVector(12);	// MagicNumber
		paramsVector(0) = this->rotationMatrix.r_11;
		paramsVector(1) = this->rotationMatrix.r_12;
		paramsVector(2) = this->rotationMatrix.r_13;
		paramsVector(3) = this->rotationMatrix.r_21;
		paramsVector(4) = this->rotationMatrix.r_22;
		paramsVector(5) = this->rotationMatrix.r_23;
		paramsVector(6) = this->f;
		paramsVector(7) = this->cameraDisplacementVector.tau_x;
		paramsVector(8) = this->cameraDisplacementVector.tau_y;
		paramsVector(9) = this->cameraDisplacementVector.tau_z;
		paramsVector(10) = this->opticalAxisVector.o_x;
		paramsVector(11) = this->opticalAxisVector.o_y;
		return paramsVector;
	}


	ProjectionParameters createFromVector(Eigen::VectorXf sourceVector)
	{
		assert(sourceVector.size() == 12);

		return ProjectionParameters(
			RotationMatrix(sourceVector(0), sourceVector(1), sourceVector(2), sourceVector(3), sourceVector(4), sourceVector(5)),
			sourceVector(6),
			CameraDisplacementVector(sourceVector(7), sourceVector(8), sourceVector(9)),
			OpticalAxisVector(sourceVector(10), sourceVector(11))
		);
	}
};