#include "box.h"
#include <cmath>

//------------------------------------------------------------------------------------------
Box::Box( const std::vector<double>& init_vec ) :
	x_min( init_vec[0] ), widht( init_vec[1] ), y_min( init_vec[2] ),
		height( init_vec[3] ), phi_min( init_vec[4] ), phi_range( init_vec[5]) {}

//------------------------------------------------------------------------------------------
void Box::GetParameters( std::vector<double>& params_vec ) const
{
	params_vec.push_back( x_min );
	params_vec.push_back( widht );
	params_vec.push_back( y_min );
	params_vec.push_back( height );
	params_vec.push_back( phi_min );
	params_vec.push_back( phi_range );
}

//------------------------------------------------------------------------------------------
void Box::GetWidhtHeight( double& par_width, double& par_height ) const
{
	par_width = widht;
	par_height = height;
}

//------------------------------------------------------------------------------------------
void Box::GetAngleRange( double& angle_min, double& angle_range ) const
{
	angle_min = phi_min;
	angle_range = phi_range;
}

//------------------------------------------------------------------------------------------
double Box::GetDiagonal() const
{
	return sqrt( widht*widht + height*height );
}