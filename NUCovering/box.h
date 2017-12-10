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
	double widht;
	/// min coordinate y
	double y_min;
	/// difference between y_min and y_max
	double height;
	/// angle 
	double phi_min;
	/// difference between phi_min and phi_max
	double phi_range;
public:
	/// default constructor
	Box() {}

	/// parametrized constructor
	Box(const std::vector<double>& init_vec);

	/// function get_box_parameters() returns all box parameters
	void GetParameters(std::vector<double>& params_vec) const;

	/// function get_widht_height() needed for getting box's widht and height
	void GetWidhtHeight(double& par_width, double& par_height) const;

	void GetAngleRange(double& angle_min, double& angle_range) const;

	/// function GetBoxDiagonal() return box's diagonal
	double GetDiagonal() const;
};