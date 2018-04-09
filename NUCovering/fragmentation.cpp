#define _USE_MATH_DEFINES
#include <cmath>
#include <engine.h>
#include <fstream>
#include <algorithm>
#include "fragmentation.h"
#include "advisor-annotate.h"

using std::vector;

//******************************************************************************************
// // ����������� ��������
//******************************************************************************************
// 2 �� ������
const double TWO_PI = M_PI * 2.0;

const double P_MIN2 = p_min * p_min;
const double P_MAX2 = p_max * p_max;

// �����������, ������������ ��� �������� �������� ���� �� �������� � �������
const double DEGREES_2_RADIAN = M_PI / 180.0;
// �����������, ������������ ��� �������� �������� ���� �� �������� � �������
const double RADIAN_2_DEGREES = 180.0 / M_PI;

// �������� �� ������ � ��������
const double PI_DEGREES = 180.0;
// �������� 2 �� ������ � ��������
const double TWO_PI_DEGREES = 360.0;

// ���������� ����� � �������� (��������� �� ���� ����)
const unsigned int SQUARED_POINTS_PER_AXIS = POINTS_PER_AXIS*POINTS_PER_AXIS;
// ����� ���������� ����� � 3-� ���������� � ������ ���������������� box'a 
const unsigned int TOTAL_POINTS_IN_BOX = SQUARED_POINTS_PER_AXIS*POINTS_PER_AXIS;

// ����� ���������� ������� gj
const unsigned int AMOUNT_OF_Gjs = 18;
//******************************************************************************************



//------------------------------------------------------------------------------------------
inline double Xci(const unsigned int& i, const double& phi)
{
	return xa[i] - xb[i] * cos(phi * DEGREES_2_RADIAN) + yb[i] * sin(phi * DEGREES_2_RADIAN);
}

//------------------------------------------------------------------------------------------
inline double Yci(const unsigned int& i, const double& phi)
{
	return ya[i] - xb[i] * sin(phi * DEGREES_2_RADIAN) - yb[i] * cos(phi * DEGREES_2_RADIAN);
}

//------------------------------------------------------------------------------------------
// ������� ���������� x ������������ �����, ������� ��
// ����������� ��������� O'xy, ������������ ����������� ��������� Oxy
inline double XbiOxy(const unsigned int& i, const double& x, const double& phi)
{
	return x + xb[i] * cos(phi * DEGREES_2_RADIAN) - yb[i] * sin(phi * DEGREES_2_RADIAN);
}

//------------------------------------------------------------------------------------------
// ������� ���������� y ������������ �����, ������� ��
// ����������� ��������� O'xy, ������������ ����������� ��������� Oxy
inline double YbiOxy(const unsigned int& i, const double& y, const double& phi)
{
	return y + xb[i] * sin(phi * DEGREES_2_RADIAN) + yb[i] * cos(phi * DEGREES_2_RADIAN);
}


//******************************************************************************************
// // ����������� ������� gj()
//******************************************************************************************
void g2i_1(const double& x, const double& y, const double& phi, double* func_values)
{
	double a1 = XbiOxy(0, x, phi) - xa[0];
	double a2 = YbiOxy(0, y, phi) - ya[0];
    func_values[0] = (a1*a1 + a2*a2 - P_MAX2);
	
	//
	a1 = XbiOxy(1, x, phi) - xa[1];
	a2 = YbiOxy(1, y, phi) - ya[1];
    func_values[1] = (a1*a1 + a2*a2 - P_MAX2);

	//
	a1 = XbiOxy(2, x, phi) - xa[2];
	a2 = YbiOxy(2, y, phi) - ya[2];
    func_values[2] = (a1*a1 + a2*a2 - P_MAX2);
}

//------------------------------------------------------------------------------------------
void g2i(const double& x, const double& y, const double& phi, double* func_values)
{
	double a1 = XbiOxy(0, x, phi) - Xci(0, phi);
	double a2 = YbiOxy(0, y, phi) - Yci(0, phi);
	func_values[3] = (P_MIN2 - a1*a1 - a2*a2);

	//
	a1 = XbiOxy(1, x, phi) - Xci(1, phi);
	a2 = YbiOxy(1, y, phi) - Yci(1, phi);
    func_values[4] = (P_MIN2 - a1*a1 - a2*a2);

	//
	a1 = XbiOxy(2, x, phi) - Xci(2, phi);
	a2 = YbiOxy(2, y, phi) - Yci(2, phi);
    func_values[5] = (P_MIN2 - a1*a1 - a2*a2);
}

//------------------------------------------------------------------------------------------
void g2i_5(const double& x, const double& y, const double& phi, double* func_values)
{
	double dist_y = YbiOxy(0, y, phi) - ya[0];
	double dist_x = XbiOxy(0, x, phi) - xa[0];
	//func_values.push_back((atan2(dist_y, dist_x) * RADIAN_2_DEGREES) - theta_ai_max[0]);
	
	////double comp = (atan2(diff_y, diff_x) * RADIAN_2_DEGREES) - theta_ai_max[0];

	// modified 27.02.18
	double rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
        func_values[6] = (acos(dist_x / rho_i) * RADIAN_2_DEGREES - theta_ai_max[0]);
	else 
        func_values[6] = ((TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES - theta_ai_max[0]);

	//
	dist_y = YbiOxy(1, y, phi) - ya[1];
	dist_x = XbiOxy(1, x, phi) - xa[1];
	//func_values.push_back((atan2(dist_y, dist_x) * RADIAN_2_DEGREES) - theta_ai_max[1]);

	////comp = (atan2(diff_y, diff_x) * RADIAN_2_DEGREES) - theta_ai_max[1];

	// modified 27.02.18
	rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
        func_values[7] = (acos(dist_x / rho_i) * RADIAN_2_DEGREES - theta_ai_max[1]);
	else 
        func_values[7] = ((TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES - theta_ai_max[1]);

	//
	dist_y = YbiOxy(2, y, phi) - ya[2];
	dist_x = XbiOxy(2, x, phi) - xa[2];
	//func_values.push_back((atan2(dist_y, dist_x) * RADIAN_2_DEGREES) - theta_ai_max[2]);

	////comp = (atan2(diff_y, diff_x) * RADIAN_2_DEGREES) - theta_ai_max[2];

	// modified 27.02.18
	rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
        func_values[8] = (acos(dist_x / rho_i) * RADIAN_2_DEGREES - theta_ai_max[2]);
	else 
        func_values[8] = ((TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES - theta_ai_max[2]);
}

//------------------------------------------------------------------------------------------
void g2i_6(const double& x, const double& y, const double& phi, double* func_values)
{
	double dist_y = YbiOxy(0, y, phi) - ya[0];
	double dist_x = XbiOxy(0, x, phi) - xa[0];
	//func_values.push_back(theta_ai_min[0] - (atan2(dist_y, dist_x) * RADIAN_2_DEGREES));

	// modified 27.02.18
	double rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);
	if (dist_y >= 0.0) 
        func_values[9] = (theta_ai_min[0] - acos(dist_x / rho_i) * RADIAN_2_DEGREES);
	else 
        func_values[9] = (theta_ai_min[0] + (acos(dist_x / rho_i) - TWO_PI) * RADIAN_2_DEGREES);

	//
	dist_y = YbiOxy(1, y, phi) - ya[1];
	dist_x = XbiOxy(1, x, phi) - xa[1];
	//func_values.push_back(theta_ai_min[1] - (atan2(dist_y, dist_x) * RADIAN_2_DEGREES));

	// modified 27.02.18
	rho_i = sqrt(dist_x *  dist_x + dist_y * dist_y);
	if (dist_y >= 0.0) 
        func_values[10] = (theta_ai_min[1] - acos(dist_x / rho_i) * RADIAN_2_DEGREES);
	else 
        func_values[10] = (theta_ai_min[1] + (acos(dist_x / rho_i) - TWO_PI) * RADIAN_2_DEGREES);

	//
	dist_y = YbiOxy(2, y, phi) - ya[2];
	dist_x = XbiOxy(2, x, phi) - xa[2];
	//func_values.push_back(theta_ai_min[2] - (atan2(dist_y, dist_x) * RADIAN_2_DEGREES));

	// modified 27.02.18
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);
	if (dist_y >= 0.0) 
        func_values[11] = (theta_ai_min[2] - acos(dist_x / rho_i) * RADIAN_2_DEGREES);
	else 
        func_values[11] = (theta_ai_min[2] + (acos(dist_x / rho_i) - TWO_PI) * RADIAN_2_DEGREES);
}

//------------------------------------------------------------------------------------------
void g2i_11(const double& x, const double& y, const double& phi, double* func_values)
{
	double dist_x = XbiOxy(0, x, phi) - xa[0];
	double dist_y = YbiOxy(0, y, phi) - ya[0];
	double rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	double theta_ai;

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
        func_values[12] = ( theta_ai - phi + PI_DEGREES - theta_bi_max[0]);
	else
        func_values[12] = (theta_ai - phi - PI_DEGREES - theta_bi_max[0]);

	//
	dist_x = XbiOxy(1, x, phi) - xa[1];
	dist_y = YbiOxy(1, y, phi) - ya[1];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
        func_values[13] = (theta_ai - phi + PI_DEGREES - theta_bi_max[1]);
	else
        func_values[13] = (theta_ai - phi - PI_DEGREES - theta_bi_max[1]);

	//
	dist_x = XbiOxy(2, x, phi) - xa[2];
	dist_y = YbiOxy(2, y, phi) - ya[2];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
        func_values[14] = (theta_ai - phi + PI_DEGREES - theta_bi_max[2]);
	else
        func_values[14] = (theta_ai - phi - PI_DEGREES - theta_bi_max[2]);
}

//------------------------------------------------------------------------------------------
void g2i_12(const double& x, const double& y, const double& phi, double* func_values)
{
	double dist_x = XbiOxy(0, x, phi) - xa[0];
	double dist_y = YbiOxy(0, y, phi) - ya[0];
	double rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	double theta_ai;

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
        func_values[15] = (theta_bi_min[0] - theta_ai + phi - PI_DEGREES);
	else
        func_values[15] = (theta_bi_min[0] - theta_ai + phi + PI_DEGREES);

	//
	dist_x = XbiOxy(1, x, phi) - xa[1];
	dist_y = YbiOxy(1, y, phi) - ya[1];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
        func_values[16] = (theta_bi_min[1] - theta_ai + phi - PI_DEGREES);
	else
        func_values[16] = (theta_bi_min[1] - theta_ai + phi + PI_DEGREES);

	//
	dist_x = XbiOxy(2, x, phi) - xa[2];
	dist_y = YbiOxy(2, y, phi) - ya[2];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
        func_values[17] = (theta_bi_min[2] - theta_ai + phi - PI_DEGREES);
	else
        func_values[17] = (theta_bi_min[2] - theta_ai + phi + PI_DEGREES);
}
//******************************************************************************************



//------------------------------------------------------------------------------------------
void GetFuncValuesInPoint(const double& x, const double& y, const double& phi, double* func_values)
{
	g2i_1(x, y, phi, func_values);
	g2i(x, y, phi, func_values);
	g2i_5(x, y, phi, func_values);
	g2i_6(x, y, phi, func_values);
	g2i_11(x, y, phi, func_values);
	g2i_12(x, y, phi, func_values);
}



/// ��������� box-��, �������� � ������������� �������� �����������
//cilk::reducer< cilk::op_vector<Box> > solution;
vector<Box> solution;
/// ��������� box-��, �� �������� � ������������� �������� �����������
//cilk::reducer< cilk::op_vector<Box> > not_solution;
vector<Box> not_solution;
/// ��������� ��������� box-��, ������� �� ��������� �� � solution, �� � not_solution
//cilk::reducer< cilk::op_vector<Box> > boundary;
vector<Box> boundary;
/// ��������� box-��, �������������� �� ��������� �������� ���������
//cilk::reducer< cilk::op_vector<Box> > next_iter_boxes;
vector<Box> next_iter_boxes;


//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(const double* init_params)
{
	current_box = Box(init_params);
}

//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(const Box& box)
{
	current_box = box;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::SplitByX(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double* box_params = new double[6];
	box.GetParameters(box_params);

    box_params[1] *= 0.5;
	new_pair_of_boxes.first = Box(box_params);

    box_params[0] += box_params[1];
	new_pair_of_boxes.second = Box(box_params);
    delete[] box_params;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::SplitByY(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double* box_params = new double[6];
	box.GetParameters(box_params);

    box_params[3] *= 0.5;
	new_pair_of_boxes.first = Box(box_params);

    box_params[2] += box_params[3];
	new_pair_of_boxes.second = Box(box_params);
    delete[] box_params;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::SplitByPhi(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double* box_params = new double[6];
	box.GetParameters(box_params);

    box_params[5] *= 0.5;
	new_pair_of_boxes.first = Box(box_params);

    box_params[4] += box_params[5];
	new_pair_of_boxes.second = Box(box_params);
    delete[] box_params;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double x_range, y_range, phi_range;
	box.GetXYPhiRanges(x_range, y_range, phi_range);

	if (x_range > y_range)
	{
		if (x_range > phi_range) SplitByX(box, new_pair_of_boxes);
		else SplitByPhi(box, new_pair_of_boxes);
	}
	else
	{
		if (y_range > phi_range) SplitByY(box, new_pair_of_boxes);
		else SplitByPhi(box, new_pair_of_boxes);
	}
}

//------------------------------------------------------------------------------------------
unsigned int low_level_fragmentation::FindTreeDepth()
{
	//puts("FindTreeDepth()");
	double box_diagonal = current_box.GetDiagonal();

	if (box_diagonal <= g_precision)
	{
		return 0;
	}
	else
	{
		boxes_pair new_boxes;
		GetNewBoxes(current_box, new_boxes);
		unsigned int tree_depth = 1;

		box_diagonal = new_boxes.first.GetDiagonal();

		if (box_diagonal <= g_precision)
		{
			return tree_depth;
		}
		else
		{
			for (;;)
			{
				GetNewBoxes(new_boxes.first, new_boxes);
				++tree_depth;

				box_diagonal = new_boxes.first.GetDiagonal();

				if (box_diagonal <= g_precision)
					break;
			}
			return tree_depth;
		}
	}
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::ClasifyBox(const Box& box, const double* func_values, unsigned& sz) const
{
	//puts("ClasifyBox()");

	if (sz == TOTAL_POINTS_IN_BOX)
	{
		// ���� ����������; ������ �������� - ��������� �� ���������� �������,
		// ������ - �� ����������.
		auto min_max_elem = std::minmax_element(&func_values[0], (func_values + TOTAL_POINTS_IN_BOX));

        if ( *min_max_elem.second < 0.0 )				// min_max_elem.second - �������� �� ��������
        {
            ANNOTATE_LOCK_ACQUIRE( 0 );
            solution.emplace_back(box);				// ��������������� box ������ �� ��������� �������
            ANNOTATE_LOCK_RELEASE( 0 );
        }
        else if ( *min_max_elem.first > 0.0 )			// min_max_elem.first - �������� �� �������
        {
            ANNOTATE_LOCK_ACQUIRE( 0 );
            not_solution.emplace_back(box);			// ��������������� box �� ������ �� ��������� �������
            ANNOTATE_LOCK_RELEASE( 0 );
        }
		else										// ��������������� box �������� ���������� � ����������� �������
		{
			boxes_pair another_boxes;
            GetNewBoxes( box, another_boxes );
            ANNOTATE_LOCK_ACQUIRE( 0 );
			next_iter_boxes.emplace_back(another_boxes.first);
			next_iter_boxes.emplace_back(another_boxes.second);
            ANNOTATE_LOCK_RELEASE( 0 );
		}	
	}
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetBoxType( const Box& box)
{
	//puts("GetBoxType()");

	// ������ �������� �������� ���������� ���� ������� gj � �������� box'a 
	double* func_values = new double[TOTAL_POINTS_IN_BOX];
    unsigned sz = 0;

	GetFusncValues(box, func_values, sz);

	ClasifyBox(box, func_values, sz);
}

//------------------------------------------------------------------------------------------
void GetGjValues(const Box& box, double* gj_values, unsigned& sz)
{
    double* box_params = new double[6];
    box.GetParameters(box_params);

    // �������� ������� �.�. ���-�� ��������, �� ������� ������� ��� ���������
    // ������ �� ������� ������, ��� ���-�� ����� �� ���� ���
    unsigned parts_per_axis = POINTS_PER_AXIS - 1;

	// ����������� ������ �� ��� x
	double x_min = box_params[0];
	// ����������� ���� ����� �� ���������� x
	double step_x = box_params[1] / parts_per_axis;

	// ����������� ������ �� ��� y
	double y_min = box_params[2];
	// ����������� ���� ����� �� ���������� y
	double step_y = box_params[3] / parts_per_axis;

	// ����������� ������ �� ��� phi
	double phi_min = box_params[4];
	// ����������� ���� ����� �� ���������� phi
	double step_phi = box_params[5] / parts_per_axis;

    delete[] box_params;

    // ������ ���������� �������� ���� 18 ������� � ������ ���������� ����� (x, y, phi)
	double* func_vals_in_curr_point = new double[AMOUNT_OF_Gjs];

    // ���������� ������ ��������������� ����� �� �����
    double x, y, phi;
    unsigned k;

    ANNOTATE_SITE_BEGIN( grid )
        // ���� �� ���� ������, ������� ������ box-�
	for (k = 0; k < TOTAL_POINTS_IN_BOX; ++k)
	{ 
        ANNOTATE_TASK_BEGIN(brut�_force)
		// ���������� ���������� ��������� ����� �� �����
		x = x_min + step_x * (k % POINTS_PER_AXIS);
		y = y_min + step_y * ((k / POINTS_PER_AXIS) % POINTS_PER_AXIS);
		phi = phi_min + step_phi * (k / SQUARED_POINTS_PER_AXIS);

		GetFuncValuesInPoint(x, y, phi, func_vals_in_curr_point);

		auto max_elem_iter = std::max_element(&func_vals_in_curr_point[0], (func_vals_in_curr_point + AMOUNT_OF_Gjs));
        gj_values[k] = *max_elem_iter;

        ANNOTATE_TASK_END(brut�_force)
	}
    ANNOTATE_SITE_END( grid )
    sz = TOTAL_POINTS_IN_BOX;

    delete[] func_vals_in_curr_point;
}

//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis(const double* init_params) :
	low_level_fragmentation(init_params) {}

//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis(Box& box) : low_level_fragmentation(box) {}

//------------------------------------------------------------------------------------------ 
void high_level_analysis::GetFusncValues(const Box& box, double* max_vals, unsigned& sz)
{
	//puts("GetMaxs()");

	if (box.GetDiagonal() > g_precision)
	{
		// ���������� ���������� ���� 18 ������� gj() �� ��������������� box-�
		GetGjValues(box, max_vals, sz);
	}
	else	// ���������� �������� �������� ������������� �������� ������������
	{
        ANNOTATE_LOCK_ACQUIRE( 0 );
		boundary.emplace_back(box);
        ANNOTATE_LOCK_RELEASE( 0 );
	}
}

//------------------------------------------------------------------------------------------
void high_level_analysis::GetSolution()
{
	unsigned int tree_depth = FindTreeDepth();

	//printf("Tree depth is %d \n", tree_depth);

    vector<Box> curr_iter_boxes;
    curr_iter_boxes.emplace_back(current_box);
    size_t number_of_box_on_level;

	puts("Getting solution");

	for (unsigned int level = 0; level < (tree_depth + 1); ++level)
	{
        number_of_box_on_level = curr_iter_boxes.size();

        ANNOTATE_SITE_BEGIN(boxes)
		for (size_t i = 0; i < number_of_box_on_level; ++i)
		{
            ANNOTATE_TASK_BEGIN(get_box_type)
			//printf("i is %d\n", i);
			GetBoxType(curr_iter_boxes[i]);
            ANNOTATE_TASK_END(get_box_type)
		}
        ANNOTATE_SITE_END(boxes)

        curr_iter_boxes = next_iter_boxes;
        // ! ���������� ������� ���������� next_iter_boxes ����� ��������� ���������
        next_iter_boxes.clear();
	}
}

void FindRingCenters(const double* cmp_angles, const unsigned& cmp_angles_size, const unsigned& rings_count, double* p_u, double* q_u)
{
	for (unsigned i = 0; i < cmp_angles_size; ++i)
	{
		double ksi = cmp_angles[i];
		double cos_ksi = cos(ksi * DEGREES_2_RADIAN);
		double sin_ksi = sin(ksi * DEGREES_2_RADIAN);

		for (unsigned j = 0; j < rings_count; ++j)
		{
			unsigned id = rings_count * i + j;

			p_u[id] = xa[j] - xb[j] * cos_ksi + yb[j] * sin_ksi;
			q_u[id] = ya[j] - xb[j] * sin_ksi - yb[j] * cos_ksi;
		}
	}
}


//------------------------------------------------------------------------------------------
void PrintWorkspace(const double* cmp_angles,const unsigned& cmp_angles_size, const unsigned& rings_count, double& phi_min, double& phi_max)
{

	const size_t solution_sz = solution.size();
	const size_t not_solution_sz = not_solution.size();
	const size_t boundary_sz = boundary.size();

	const unsigned number_of_params = 6;

	printf("Solution sz = %Iu\n", solution_sz);
	printf("NOT Solution sz = %Iu\n", not_solution_sz);
	printf("Boundary sz = %Iu\n", boundary_sz);

	double params[6];

	solution[0].GetParameters(params);
	printf("\nSol[0]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);
    solution[1].GetParameters(params);
    printf("\nSol[1]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);
    solution[2].GetParameters(params);
    printf("\nSol[2]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);

	not_solution[0].GetParameters(params);
	printf("Not_sol[0]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);
    not_solution[1].GetParameters(params);
    printf("Not_sol[1]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);
    not_solution[2].GetParameters(params);
    printf("Not_sol[2]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);

	boundary[0].GetParameters(params);
	printf("Boundary[0]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);
    boundary[1].GetParameters(params);
    printf("Boundary[1]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);
    boundary[2].GetParameters(params);
    printf("Boundary[2]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);

	// solution
	double* sol_arr = new double[number_of_params*solution_sz];
	vector<Box>::const_iterator beg_iter = solution.begin();
	vector<Box>::const_iterator end_iter = solution.end();

	unsigned i = 0;
	for (beg_iter; beg_iter != end_iter; ++beg_iter, ++i)
	{
		beg_iter->GetParameters(params);
		memcpy(sol_arr + i*number_of_params, params, number_of_params * sizeof(double));
	}

	// not_soltion
	double* not_sol_arr = new double[number_of_params*not_solution_sz];
	beg_iter = not_solution.begin();
	end_iter = not_solution.end();

	i = 0;
	for (beg_iter; beg_iter != end_iter; ++beg_iter, ++i)
	{
		beg_iter->GetParameters(params);
		memcpy(not_sol_arr + i*number_of_params, params, number_of_params * sizeof(double));
	}

	// boundary
	double* boundary_arr = new double[number_of_params*boundary_sz];
	beg_iter = boundary.begin();
	end_iter = boundary.end();

	i = 0;
	for (beg_iter; beg_iter != end_iter; ++beg_iter, ++i)
	{
		beg_iter->GetParameters(params);
		memcpy(boundary_arr + i*number_of_params, params, number_of_params * sizeof(double));
	}

	// ������������� ������ � VS
	solution.shrink_to_fit();
	not_solution.shrink_to_fit();
	boundary.shrink_to_fit();

	// ������ ������  Matlab 
	puts("Printing workspace!");

	double* p_u = new double[rings_count * cmp_angles_size];
	double* q_u = new double[rings_count * cmp_angles_size];

	FindRingCenters(cmp_angles, cmp_angles_size, rings_count, p_u, q_u);

	double vs_R[] = { p_min, p_max };
	double vs_phi_min_max[] = { phi_min, phi_max };

	// ������ ������ MATLAB
	Engine* engine_ptr = engOpen(NULL);

	mxArray* matlab_xc				= mxCreateDoubleMatrix(rings_count, cmp_angles_size, mxREAL);
	mxArray* matlab_yc				= mxCreateDoubleMatrix(rings_count, cmp_angles_size, mxREAL);
	mxArray* matlab_R				= mxCreateDoubleMatrix(1, 2, mxREAL);
	mxArray* matlab_phi_min_max		= mxCreateDoubleMatrix(1, 2, mxREAL);

	mxArray* matlab_sol_array		= mxCreateDoubleMatrix(number_of_params, solution_sz, mxREAL);
	mxArray* matlab_not_sol_array	= mxCreateDoubleMatrix(number_of_params, not_solution_sz, mxREAL);
	mxArray* matlab_boundary_array	= mxCreateDoubleMatrix(number_of_params, boundary_sz, mxREAL);

	mxArray* matlab_cmp_angles		= mxCreateDoubleMatrix(1, cmp_angles_size, mxREAL);
	mxArray* matlab_delta			= mxCreateDoubleMatrix(1, 1, mxREAL);

	memcpy(mxGetPr(matlab_delta), &g_precision, sizeof(double));
	engPutVariable(engine_ptr, "delta", matlab_delta);

	memcpy(mxGetPr(matlab_cmp_angles), &cmp_angles[0], cmp_angles_size * sizeof(double));
	engPutVariable(engine_ptr, "cmp_angles", matlab_cmp_angles);

	memcpy(mxGetPr(matlab_xc), &p_u[0], rings_count * cmp_angles_size * sizeof(double));
	engPutVariable(engine_ptr, "xc", matlab_xc);
	engEvalString(engine_ptr, "xc = xc';");				// ���������� ������������� ���������������!

	memcpy(mxGetPr(matlab_yc), &q_u[0], rings_count * cmp_angles_size * sizeof(double));
	engPutVariable(engine_ptr, "yc", matlab_yc);
	engEvalString(engine_ptr, "yc = yc';");				// ���������� ������������� ���������������!

	memcpy(mxGetPr(matlab_R), &vs_R[0], 2 * sizeof(double));
	engPutVariable(engine_ptr, "R", matlab_R);

	memcpy(mxGetPr(matlab_phi_min_max), &vs_phi_min_max[0], 2 * sizeof(double));
	engPutVariable(engine_ptr, "phi_min_max", matlab_phi_min_max);

	// �������� ������ box'��, ���������� �������� �������� ������� 
	memcpy(mxGetPr(matlab_sol_array), &sol_arr[0], number_of_params * solution_sz * sizeof(double));
	engPutVariable(engine_ptr, "sol_array", matlab_sol_array);
	engEvalString(engine_ptr, "sol_array = sol_array';");				// ���������� ������������� ���������������!

	// �������� ������ box'��, �� ���������� �������� �������� �������
	memcpy(mxGetPr(matlab_not_sol_array), &not_sol_arr[0], number_of_params * not_solution_sz * sizeof(double));
	engPutVariable(engine_ptr, "not_sol_array", matlab_not_sol_array);
	engEvalString(engine_ptr, "not_sol_array = not_sol_array';");		// ���������� ������������� ���������������!

	// �������� ������ box'��, ������� �� �������
	memcpy(mxGetPr(matlab_boundary_array), &boundary_arr[0], number_of_params * boundary_sz * sizeof(double));
	engPutVariable(engine_ptr, "boundary", matlab_boundary_array);
	engEvalString(engine_ptr, "boundary = boundary';");					// ���������� ������������� ���������������!

	// ������������� ������ � VS
	delete[] sol_arr;
	delete[] not_sol_arr;
	delete[] boundary_arr;

	const char* plotting_script_dir = "cd D:\\Study\\Master\\Plotting";
	engEvalString(engine_ptr, plotting_script_dir);
	engEvalString(engine_ptr, "PrintWorkspace(xc, yc, R, phi_min_max, sol_array, not_sol_array, boundary, cmp_angles, delta)");
	 
	// ������������� ������
	//mxDestroyArray(matlab_xc);
	//mxDestroyArray(matlab_yc);
	//mxDestroyArray(matlab_R);
	//mxDestroyArray(matlab_phi_min_max);
	//mxDestroyArray(matlab_sol_array);
	//mxDestroyArray(matlab_not_sol_array);
	//mxDestroyArray(matlab_boundary_array);
	//mxDestroyArray(matlab_cmp_angles);
	//mxDestroyArray(matlab_delta);

	//engClose(engine_ptr);
}
