#define _USE_MATH_DEFINES
#include <cmath>
#include <Engine.h>
//#include <cilk/cilk.h>
//#include <cilk/reducer_vector.h>
#include <fstream>
#include <algorithm>
#include "fragmentation.h"

using std::vector;

//******************************************************************************************
// // Определение констант
//******************************************************************************************
// 2 ПИ радиан
const double TWO_PI = M_PI * 2.0;

const double P_MIN2 = p_min*p_min;
const double P_MAX2 = p_max*p_max;

// коэффициент, используемый для перевода значения угла из градусов в радианы
const double DEGREES_2_RADIAN = M_PI / 180.0;
// коэффициент, используемый для перевода значения угла из радианов в градусы
const double RADIAN_2_DEGREES = 180.0 / M_PI;

// значение ПИ радиан в градусах
const double PI_DEGREES = 180.0;
// значение 2 ПИ радиан в градусах
const double TWO_PI_DEGREES = 360.0;

// количество точек в квадрате (суммарное по двум осям)
const unsigned int SQUARED_POINTS_PER_AXIS = points_per_axis*points_per_axis;
// общее количество точек в 3-х измерениях в рамках рассматриваемого box'a 
const unsigned int TOTAL_POINTS_IN_BOX = SQUARED_POINTS_PER_AXIS*points_per_axis;

// общее количество функций gj
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
// находим координату x произвольной точки, лежащей во
// вращающейся плоскости O'xy, относительно неподвижной плоскости Oxy
inline double XbiOxy(const unsigned int& i, const double& x, const double& phi)
{
	return x + xb[i] * cos(phi * DEGREES_2_RADIAN) - yb[i] * sin(phi * DEGREES_2_RADIAN);
}

//------------------------------------------------------------------------------------------
// находим координату y произвольной точки, лежащей во
// вращающейся плоскости O'xy, относительно неподвижной плоскости Oxy
inline double YbiOxy(const unsigned int& i, const double& y, const double& phi)
{
	return y + xb[i] * sin(phi * DEGREES_2_RADIAN) + yb[i] * cos(phi * DEGREES_2_RADIAN);
}


//******************************************************************************************
// // Определение функций gj()
//******************************************************************************************
void g2i_1(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double a1 = XbiOxy(0, x, phi) - xa[0];
	double a2 = YbiOxy(0, y, phi) - ya[0];
	func_values.push_back(a1*a1 + a2*a2 - P_MAX2);
	
	//
	a1 = XbiOxy(1, x, phi) - xa[1];
	a2 = YbiOxy(1, y, phi) - ya[1];
	func_values.push_back(a1*a1 + a2*a2 - P_MAX2);

	//
	a1 = XbiOxy(2, x, phi) - xa[2];
	a2 = YbiOxy(2, y, phi) - ya[2];
	func_values.push_back(a1*a1 + a2*a2 - P_MAX2);
}

//------------------------------------------------------------------------------------------
void g2i(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double a1 = XbiOxy(0, x, phi) - Xci(0, phi);
	double a2 = YbiOxy(0, y, phi) - Yci(0, phi);
	func_values.push_back(P_MIN2 - a1*a1 - a2*a2);

	//
	a1 = XbiOxy(1, x, phi) - Xci(1, phi);
	a2 = YbiOxy(1, y, phi) - Yci(1, phi);
	func_values.push_back(P_MIN2 - a1*a1 - a2*a2);

	//
	a1 = XbiOxy(2, x, phi) - Xci(2, phi);
	a2 = YbiOxy(2, y, phi) - Yci(2, phi);
	func_values.push_back(P_MIN2 - a1*a1 - a2*a2);
}

//------------------------------------------------------------------------------------------
void g2i_5(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double dist_y = YbiOxy(0, y, phi) - ya[0];
	double dist_x = XbiOxy(0, x, phi) - xa[0];
	//func_values.push_back((atan2(dist_y, dist_x) * RADIAN_2_DEGREES) - theta_ai_max[0]);
	
	////double comp = (atan2(diff_y, diff_x) * RADIAN_2_DEGREES) - theta_ai_max[0];

	// modified 27.02.18
	double rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
		func_values.push_back(acos(dist_x / rho_i)*RADIAN_2_DEGREES - theta_ai_max[0]);
	else 
		func_values.push_back((TWO_PI - acos(dist_x / rho_i))*RADIAN_2_DEGREES - theta_ai_max[0]);

	//
	dist_y = YbiOxy(1, y, phi) - ya[1];
	dist_x = XbiOxy(1, x, phi) - xa[1];
	//func_values.push_back((atan2(dist_y, dist_x) * RADIAN_2_DEGREES) - theta_ai_max[1]);

	////comp = (atan2(diff_y, diff_x) * RADIAN_2_DEGREES) - theta_ai_max[1];

	// modified 27.02.18
	rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
		func_values.push_back(acos(dist_x / rho_i)*RADIAN_2_DEGREES - theta_ai_max[1]);
	else 
		func_values.push_back((TWO_PI - acos(dist_x / rho_i))*RADIAN_2_DEGREES - theta_ai_max[1]);

	//
	dist_y = YbiOxy(2, y, phi) - ya[2];
	dist_x = XbiOxy(2, x, phi) - xa[2];
	//func_values.push_back((atan2(dist_y, dist_x) * RADIAN_2_DEGREES) - theta_ai_max[2]);

	////comp = (atan2(diff_y, diff_x) * RADIAN_2_DEGREES) - theta_ai_max[2];

	// modified 27.02.18
	rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
		func_values.push_back(acos(dist_x / rho_i)*RADIAN_2_DEGREES - theta_ai_max[2]);
	else 
		func_values.push_back((TWO_PI - acos(dist_x / rho_i))*RADIAN_2_DEGREES - theta_ai_max[2]);
}

//------------------------------------------------------------------------------------------
void g2i_6(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double dist_y = YbiOxy(0, y, phi) - ya[0];
	double dist_x = XbiOxy(0, x, phi) - xa[0];
	//func_values.push_back(theta_ai_min[0] - (atan2(dist_y, dist_x) * RADIAN_2_DEGREES));

	// modified 27.02.18
	double rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
		func_values.push_back(theta_ai_min[0] - acos(dist_x / rho_i)*RADIAN_2_DEGREES);
	else 
		func_values.push_back(theta_ai_min[0] + (acos(dist_x / rho_i) - TWO_PI)*RADIAN_2_DEGREES);

	//
	dist_y = YbiOxy(1, y, phi) - ya[1];
	dist_x = XbiOxy(1, x, phi) - xa[1];
	//func_values.push_back(theta_ai_min[1] - (atan2(dist_y, dist_x) * RADIAN_2_DEGREES));

	// modified 27.02.18
	rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
		func_values.push_back(theta_ai_min[1] - acos(dist_x / rho_i)*RADIAN_2_DEGREES);
	else 
		func_values.push_back(theta_ai_min[1] + (acos(dist_x / rho_i) - TWO_PI)*RADIAN_2_DEGREES);

	//
	dist_y = YbiOxy(2, y, phi) - ya[2];
	dist_x = XbiOxy(2, x, phi) - xa[2];
	//func_values.push_back(theta_ai_min[2] - (atan2(dist_y, dist_x) * RADIAN_2_DEGREES));

	// modified 27.02.18
	rho_i = sqrt(dist_x*dist_x + dist_y*dist_y);
	if (dist_y >= 0.0) 
		func_values.push_back(theta_ai_min[2] - acos(dist_x / rho_i)*RADIAN_2_DEGREES);
	else 
		func_values.push_back(theta_ai_min[2] + (acos(dist_x / rho_i) - TWO_PI)*RADIAN_2_DEGREES);
}

//------------------------------------------------------------------------------------------
void g2i_11(const double& x, const double& y, const double& phi, vector<double>& func_values)
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
		func_values.push_back( theta_ai - phi + PI_DEGREES - theta_bi_max[0]);
	else
		func_values.push_back(theta_ai - phi - PI_DEGREES - theta_bi_max[0]);

	//
	dist_x = XbiOxy(1, x, phi) - xa[1];
	dist_y = YbiOxy(1, y, phi) - ya[1];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
		func_values.push_back(theta_ai - phi + PI_DEGREES - theta_bi_max[1]);
	else
		func_values.push_back(theta_ai - phi - PI_DEGREES - theta_bi_max[1]);

	//
	dist_x = XbiOxy(2, x, phi) - xa[2];
	dist_y = YbiOxy(2, y, phi) - ya[2];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
		func_values.push_back(theta_ai - phi + PI_DEGREES - theta_bi_max[2]);
	else
		func_values.push_back(theta_ai - phi - PI_DEGREES - theta_bi_max[2]);
}

//------------------------------------------------------------------------------------------
void g2i_12(const double& x, const double& y, const double& phi, vector<double>& func_values)
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
		func_values.push_back(theta_bi_min[0] - theta_ai + phi - PI_DEGREES);
	else
		func_values.push_back(theta_bi_min[0] - theta_ai + phi + PI_DEGREES);

	//
	dist_x = XbiOxy(1, x, phi) - xa[1];
	dist_y = YbiOxy(1, y, phi) - ya[1];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
		func_values.push_back(theta_bi_min[1] - theta_ai + phi - PI_DEGREES);
	else
		func_values.push_back(theta_bi_min[1] - theta_ai + phi + PI_DEGREES);

	//
	dist_x = XbiOxy(2, x, phi) - xa[2];
	dist_y = YbiOxy(2, y, phi) - ya[2];
	rho_i = sqrt(dist_x * dist_x + dist_y * dist_y);

	if (dist_y >= 0.0)
		theta_ai = acos(dist_x / rho_i) * RADIAN_2_DEGREES;
	else
		theta_ai = (TWO_PI - acos(dist_x / rho_i)) * RADIAN_2_DEGREES;

	if ((theta_ai - phi + PI_DEGREES) < TWO_PI_DEGREES)
		func_values.push_back(theta_bi_min[2] - theta_ai + phi - PI_DEGREES);
	else
		func_values.push_back(theta_bi_min[2] - theta_ai + phi + PI_DEGREES);
}
//******************************************************************************************



//------------------------------------------------------------------------------------------
void GetFuncValuesInPoint(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	g2i_1(x, y, phi, func_values);
	g2i(x, y, phi, func_values);
	g2i_5(x, y, phi, func_values);
	g2i_6(x, y, phi, func_values);
	g2i_11(x, y, phi, func_values);
	g2i_12(x, y, phi, func_values);
}



/// контейнер box-ов, входящих в аппроксимацию рабочего простраства
//cilk::reducer< cilk::op_vector<Box> > solution;
vector<Box> solution;
/// контейнер box-ов, не входящих в аппроксимацию рабочего простраства
//cilk::reducer< cilk::op_vector<Box> > not_solution;
vector<Box> not_solution;
/// контейнер граничных box-ов, которые не относятся ни к solution, ни к not_solution
//cilk::reducer< cilk::op_vector<Box> > boundary;
vector<Box> boundary;
/// контейнер box-ов, обрабатываемых на следующей итерации алгоритма
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

	double x_min = box_params[0];
	double new_x_range = box_params[1] * 0.5;

	double new_box_params[6] = { x_min, new_x_range, box_params[2], box_params[3], box_params[4], box_params[5] };
	new_pair_of_boxes.first = Box(new_box_params);

	new_box_params[0] = x_min + new_x_range;
	new_pair_of_boxes.second = Box(new_box_params);
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::SplitByY(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double* box_params = new double[6];
	box.GetParameters(box_params);

	double y_min = box_params[2];
	double new_y_range = box_params[3] * 0.5;

	double new_box_params[6] = { box_params[0], box_params[1], y_min, new_y_range, box_params[4], box_params[5] };
	new_pair_of_boxes.first = Box(new_box_params);

	new_box_params[2] = y_min + new_y_range;
	new_pair_of_boxes.second = Box(new_box_params);
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::SplitByPhi(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double* box_params = new double[6];
	box.GetParameters(box_params);

	double phi_min = box_params[4];
	double new_phi_range = box_params[5] * 0.5;

	double new_box_params[6] = { box_params[0], box_params[1], box_params[2], box_params[3], phi_min, new_phi_range };
	new_pair_of_boxes.first = Box(new_box_params);

	new_box_params[4] = phi_min + new_phi_range;
	new_pair_of_boxes.second = Box(new_box_params);
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double x_range, y_range, phi_range;
	box.GetXYPhiRanges(x_range, y_range, phi_range);

	if (x_range >= y_range)
	{
		if (x_range >= phi_range) SplitByX(box, new_pair_of_boxes);
		else SplitByPhi(box, new_pair_of_boxes);
	}
	else
	{
		if (y_range >= phi_range) SplitByY(box, new_pair_of_boxes);
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
		//SplitByY(current_box, new_boxes);
		GetNewBoxes(current_box, new_boxes);
		unsigned int tree_depth = 1;

		box_diagonal = new_boxes.first.GetDiagonal();

		if (box_diagonal <= g_precision)
		{
			return tree_depth;
		}
		else
		{
			for (;; )
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
void low_level_fragmentation::ClasifyBox(const Box& box, const std::vector<double>& func_values) const
{
	puts("ClasifyBox()");

	if (func_values.empty() == false)
	{
		// Пара итераторов; первый итератор - указывает на наименьший элемент,
		// второй - на наибольший.
		auto min_max_elem = std::minmax_element(func_values.begin(), func_values.end());

		if (*min_max_elem.second < 0.0)				// min_max_elem.second - итератор на максимум
			solution.push_back(box);				// рассматриваемый box входит во множество решений
		else if (*min_max_elem.first > 0.0)			// min_max_elem.first - итератор на минимум
			not_solution.push_back(box);			// рассматриваемый box не входит во множество решений
		else										// рассматриваемый box подлежит разбинению и дальнейшему анализу
		{
			boxes_pair another_boxes;
			GetNewBoxes(box, another_boxes);				

			//next_iter_boxes->push_back(another_boxes.first);
			next_iter_boxes.push_back(another_boxes.first);
			//next_iter_boxes->push_back(another_boxes.second);
			next_iter_boxes.push_back(another_boxes.second);
		}	
	}
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetBoxType( const Box& box)
{
	//puts("GetBoxType()");

	// вектор хранящий значения максимумов всех функций gj в пределах box'a 
	vector<double> func_values;
	func_values.reserve(TOTAL_POINTS_IN_BOX);

	GetFusncValues(box, func_values);

	ClasifyBox(box, func_values);
}

//------------------------------------------------------------------------------------------
void GetGjValues(const double* box_params, vector<double>& gj_values)
{
	// определение границ по оси x
	double x_min = box_params[0];
	double x_range = box_params[1];
	double x_max = x_min + x_range;

	// вычитаем единицу т.к. кол-во отрезков, на которые делится ось координат 
	// всегда на единицу меньше, чем кол-во точек на этой оси
	unsigned parts_per_axis = points_per_axis - 1;

	// определение шага сетки по координате x
	double step_x = x_range / parts_per_axis;

	// определение границ по оси y
	double y_min = box_params[2];
	double y_range = box_params[3];
	double y_max = y_min + y_range;

	// определение шага сетки по координате y
	double step_y = y_range / parts_per_axis;

	// определение границ по оси phi
	double phi_min = box_params[4];
	double phi_range = box_params[5];
	double phi_max = phi_min + phi_range;

	// определение шага сетки по координате phi
	double step_phi = phi_range / parts_per_axis;

	// координаты первой рассматриваемой точки на сетке
	double x = x_min, y = y_min, phi = phi_min;
	
	vector<double> func_vals_in_curr_point;
	func_vals_in_curr_point.reserve(AMOUNT_OF_Gjs);

	vector<double>::const_iterator max_elem_iter;

	// цикл по всем точкам, лежащим внутри box-а
	for (unsigned k = 0; k < TOTAL_POINTS_IN_BOX; ++k)
	{ 
		// определяем координаты очередной точки на сетке
		x = x_min + step_x*(k % points_per_axis);
		y = y_min + step_y*((k / points_per_axis) % points_per_axis);
		phi = phi_min + step_phi*(k / SQUARED_POINTS_PER_AXIS);

		GetFuncValuesInPoint(x, y, phi, func_vals_in_curr_point);

		max_elem_iter = std::max_element(func_vals_in_curr_point.begin(), func_vals_in_curr_point.end());
		gj_values.push_back(*max_elem_iter);

		func_vals_in_curr_point.clear();
	}
}

//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis(const double* init_params) :
	low_level_fragmentation(init_params) {}

//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis(Box& box) : low_level_fragmentation(box) {}

//------------------------------------------------------------------------------------------ 
void high_level_analysis::GetFusncValues(const Box& box, vector<double>& max_vals)
{
	//puts("GetMaxs()");

	double curr_box_diagonal = box.GetDiagonal();

	if (curr_box_diagonal > g_precision)
	{
		double* box_params = new double[6];
		box.GetParameters(box_params);

		// нахождение максимумов всех 18 функций gj() на рассматриваемом box-е
		GetGjValues(box_params, max_vals);
	}
	else	// достигнута желаемая точность аппроксимации рабочего пространства
	{
		boundary.push_back(box);
	}
}

//------------------------------------------------------------------------------------------
void high_level_analysis::GetSolution()
{
	unsigned int tree_depth = FindTreeDepth();

	//printf("Tree depth is %d \n", tree_depth);

	next_iter_boxes.push_back(current_box);

	vector<Box> curr_iter_boxes, empty_vec;

	puts("Getting solution");

	for (unsigned int level = 0; level < (tree_depth + 1); ++level)
	{
		size_t number_of_box_on_level = next_iter_boxes.size();

		//next_iter_boxes.move_out(curr_iter_boxes);

		curr_iter_boxes = next_iter_boxes;

		// ! необходимо очищать содержимое next_iter_boxes
		//next_iter_boxes.set_value(empty_vec);

		next_iter_boxes.clear();

		//cilk_for(size_t i = 0; i < number_of_box_on_level; ++i)
		for (size_t i = 0; i < number_of_box_on_level; ++i)
		{
			//printf("i is %d\n", i);
			GetBoxType(curr_iter_boxes[i]);
		}
	}
}

//------------------------------------------------------------------------------------------
void PrintWorkspace(const double* cmp_angles, unsigned cmp_angles_size)
{
	const unsigned number_of_params = 6;

	size_t solution_sz = solution.size();
	size_t not_solution_sz = not_solution.size();
	size_t boundary_sz = boundary.size();

	printf("Solution sz = %Iu\n", solution_sz);
	printf("NOT Solution sz = %Iu\n", not_solution_sz);
	printf("Boundary sz = %Iu\n", boundary_sz);

	double params[6];

	//solution[0].GetParameters(params);
	//printf("\nSol[0]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);

	//not_solution[0].GetParameters(params);
	//printf("Not_sol[0]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);

	//boundary[0].GetParameters(params);
	//printf("Boundary[0]: %lf  %lf  %lf  %lf  %lf  %lf\n\n", params[0], params[1], params[2], params[3], params[4], params[5]);

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

	// высвобождение памяти в VS
	solution.shrink_to_fit();
	not_solution.shrink_to_fit();
	boundary.shrink_to_fit();

	// начало работы  Matlab 
	puts("Printing workspace!");

	Engine* engine_ptr = engOpen(NULL);

	mxArray* matlab_sol_array		= mxCreateDoubleMatrix(number_of_params, solution_sz, mxREAL);
	mxArray* matlab_not_sol_array	= mxCreateDoubleMatrix(number_of_params, not_solution_sz, mxREAL);
	mxArray* matlab_boundary_array	= mxCreateDoubleMatrix(number_of_params, boundary_sz, mxREAL);
	mxArray* matlab_cmp_angles		= mxCreateDoubleMatrix(1, cmp_angles_size, mxREAL);
	mxArray* matlab_delta			= mxCreateDoubleMatrix(1, 1, mxREAL);

	memcpy(mxGetPr(matlab_delta), &g_precision, sizeof(double));
	engPutVariable(engine_ptr, "delta", matlab_delta);

	memcpy(mxGetPr(matlab_cmp_angles), &cmp_angles[0], cmp_angles_size * sizeof(double));
	engPutVariable(engine_ptr, "cmp_angles", matlab_cmp_angles);

	// копируем вектор box'ов, являющихся решением исходной системы 
	memcpy(mxGetPr(matlab_sol_array), &sol_arr[0], number_of_params * solution_sz * sizeof(double));
	engPutVariable(engine_ptr, "sol_array", matlab_sol_array);
	engEvalString(engine_ptr, "sol_array = sol_array';");				// необходимо дополнительно транспонировать!

	// копируем вектор box'ов, НЕ являющихся решением исходной системы
	memcpy(mxGetPr(matlab_not_sol_array), &not_sol_arr[0], number_of_params * not_solution_sz * sizeof(double));
	engPutVariable(engine_ptr, "not_sol_array", matlab_not_sol_array);
	engEvalString(engine_ptr, "not_sol_array = not_sol_array';");		// необходимо дополнительно транспонировать!

	// копируем вектор box'ов, лежащих на границе
	memcpy(mxGetPr(matlab_boundary_array), &boundary_arr[0], number_of_params * boundary_sz * sizeof(double));
	engPutVariable(engine_ptr, "boundary", matlab_boundary_array);
	engEvalString(engine_ptr, "boundary = boundary';");					// необходимо дополнительно транспонировать!

	// высвобождение памяти в VS
	delete[] sol_arr;
	delete[] not_sol_arr;
	delete[] boundary_arr;

	const char* plotting_script_dir = "cd D:\\Study\\Master";
	engEvalString(engine_ptr, plotting_script_dir);
	engEvalString(engine_ptr, "PrintWorkspace(sol_array, not_sol_array, boundary, cmp_angles, delta)");
	 
	// высвобождение памяти
	mxDestroyArray(matlab_sol_array);
	mxDestroyArray(matlab_not_sol_array);
	mxDestroyArray(matlab_boundary_array);
	mxDestroyArray(matlab_cmp_angles);
	mxDestroyArray(matlab_delta);

	engClose(engine_ptr);
}
