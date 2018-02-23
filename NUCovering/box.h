#pragma once
#include <vector>

/// struct Box defines box parameters
/// and functions for manipulating with them
struct Box
{
protected:
	/// min coordinate õ
	double x_min;
	/// difference between x_max and x_min
	double x_range;
	/// min coordinate y
	double y_min;
	/// difference between y_min and y_max
	double y_range;
	/// angle 
	double phi_min;
	/// difference between phi_min and phi_max
	double phi_range;
public:
	/// default constructor
	Box() {}

	/// parametrized constructor
	Box(const double* iparams);

	/// function get_box_parameters() returns all box parameters
	void GetParameters(double* oparams) const;

	/// function get_widht_height() needed for getting box's widht and height
	void GetXYPhiRanges(double& ox_range, double& oy_range, double& ophi_range) const;

	void GetAngleRange(double& oangle_min, double& oangle_range) const;

	/// function GetBoxDiagonal() return box's diagonal
	double GetDiagonal() const;

	Box& operator=(const Box& other) = default;
};