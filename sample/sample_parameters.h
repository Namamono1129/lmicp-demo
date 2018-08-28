#pragma once
#include "projection_parameters.h"

ProjectionParameters GetSampleProjectionParameters()
{
	/*** Rotation Matrix ***/
	float r_11 = 1;
	float r_12 = 0;
	float r_13 = 0;
	float r_21 = 0;
	float r_22 = 1;
	float r_23 = 0;
	RotationMatrix R = RotationMatrix(r_11, r_12, r_13, r_21, r_22, r_23);

	/*** Focal length ***/
	float f = 100;
	ProjectionMatrix P = ProjectionMatrix(f);

	/*** Camera Displacement Vector ***/
	float tau_x = 0;
	float tau_y = 0;
	float tau_z = -50;
	CameraDisplacementVector vec_tau = CameraDisplacementVector(tau_x, tau_y, tau_z);

	/*** Optical Axis of projected plane ***/
	float o_x = 320;
	float o_y = 240;
	OpticalAxisVector vec_o = OpticalAxisVector(o_x, o_y);

	return ProjectionParameters(R, f, vec_tau, vec_o);
}
