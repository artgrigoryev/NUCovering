#include "box.h"
#include <cmath>

//------------------------------------------------------------------------------------------
Box::Box(const double* iparams) :
	x_min(iparams[0]), x_range(iparams[1]), y_min(iparams[2]),
	y_range(iparams[3]), phi_min(iparams[4]), phi_range(iparams[5]) {}

//------------------------------------------------------------------------------------------
void Box::GetParameters(double* oparams) const
{
	oparams[0]= x_min;
	oparams[1] = x_range;
	oparams[2] = y_min;
	oparams[3] = y_range;
	oparams[4] = phi_min;
	oparams[5] = phi_range;
}

//------------------------------------------------------------------------------------------
void Box::GetXYPhiRanges(double& ox_range, double& oy_range, double& ophi_range) const
{
	ox_range = x_range;
	oy_range = y_range;
	ophi_range = phi_range;
}

//------------------------------------------------------------------------------------------
void Box::GetAngleRange(double& oangle_min, double& oangle_range) const
{
	oangle_min = phi_min;
	oangle_range = phi_range;
}

//------------------------------------------------------------------------------------------
double Box::GetDiagonal() const
{
	return sqrt(x_range*x_range + y_range*y_range + phi_range*phi_range);
}

//------------------------------------------------------------------------------------------
//Box& Box::operator= (const Box& other)
//{
//	if (this == &other) 
//		return *this;
//
//	x_min	= other.x_min;
//	x_range = other.x_range;
//	y_min	= other.y_min;
//	y_range = other.y_range;
//	phi_min = other.phi_min;
//	phi_min = other.phi_range;
//
//	return *this;
//}