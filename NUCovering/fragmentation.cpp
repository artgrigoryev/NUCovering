#define _USE_MATH_DEFINES
#include <cmath>
#include <Engine.h>
//#include <cilk/cilk.h>
//#include <cilk/reducer_vector.h>
#include <fstream>
#include <algorithm>
#include "fragmentation.h"

using std::vector;

const double p_min2 = p_min*p_min;
const double p_max2 = p_max*p_max;
const double TWO_PI = M_PI * 2;

// общее количество точек на всех 3-х осях координат
const unsigned int points_in_3d = (unsigned) pow(points_per_axis, 3.0);

const unsigned int squared_points_per_axis = points_per_axis*points_per_axis;

// общее количество функций gj
const unsigned int amount_of_gjs = 18;
// количество функций gj в группе gj[]
const unsigned int func_per_group = 3;

// коэффициент, используемый для перевода значения угла из градусов в радианы
const double val_to_radian = M_PI / 180.0;
// коэффициент, используемый для перевода значения угла из радианов в градусы
const double radian_to_val = 180.0 / M_PI;

//------------------------------------------------------------------------------------------
inline double Xci(const unsigned int& i, const double& phi)
{
	return xa[i] - xb[i] * cos(phi * val_to_radian) + yb[i] * sin(phi * val_to_radian);
}

//------------------------------------------------------------------------------------------
inline double Yci(const unsigned int& i, const double& phi)
{
	return ya[i] - xb[i] * sin(phi * val_to_radian) - yb[i] * cos(phi* val_to_radian);
}

//------------------------------------------------------------------------------------------
// находим координату x произвольной точки, лежащей во
// вращающейся плоскости O'xy, относительно неподвижной плоскости Oxy
inline double XbiOxy(const unsigned int& i, const double& x, const double& phi)
{
	return x + xb[i] * cos(phi * val_to_radian) - yb[i] * sin(phi * val_to_radian);
}

//------------------------------------------------------------------------------------------
// находим координату y произвольной точки, лежащей во
// вращающейся плоскости O'xy, относительно неподвижной плоскости Oxy
inline double YbiOxy(const unsigned int& i, const double& y, const double& phi)
{
	return y + xb[i] * sin(phi* val_to_radian) + yb[i] * cos(phi * val_to_radian);
}

/// functions gj()
//------------------------------------------------------------------------------------------
void g2i_1(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double a1 = XbiOxy(0, x, y) - xa[0];
	double a2 = YbiOxy(0, y, phi) - ya[0];
	func_values.push_back(a1*a1 + a2*a2 - p_max2);
	
	a1 = XbiOxy(1, x, y) - xa[1];
	a2 = YbiOxy(1, y, phi) - ya[1];
	func_values.push_back(a1*a1 + a2*a2 - p_max2);

	a1 = XbiOxy(2, x, y) - xa[2];
	a2 = YbiOxy(2, y, phi) - ya[2];
	func_values.push_back(a1*a1 + a2*a2 - p_max2);
}

//------------------------------------------------------------------------------------------
void g2i(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double a1 = XbiOxy(0, x, phi) - Xci(0, phi);
	double a2 = YbiOxy(0, y, phi) - Yci(0, phi);
	func_values.push_back(p_min2 - a1*a1 - a2*a2);

	a1 = XbiOxy(1, x, phi) - Xci(1, phi);
	a2 = YbiOxy(1, y, phi) - Yci(1, phi);
	func_values.push_back(p_min2 - a1*a1 - a2*a2);

	a1 = XbiOxy(2, x, phi) - Xci(2, phi);
	a2 = YbiOxy(2, y, phi) - Yci(2, phi);
	func_values.push_back(p_min2 - a1*a1 - a2*a2);
}

//------------------------------------------------------------------------------------------
void g2i_5(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double a1 = YbiOxy(0, y, phi) - ya[0];
	double a2 = XbiOxy(0, x, phi) - xa[0];
	func_values.push_back((atan2(a1, a2) * radian_to_val) - theta_ai_max[0]);

	a1 = YbiOxy(1, y, phi) - ya[1];
	a2 = XbiOxy(1, x, phi) - xa[1];
	func_values.push_back((atan2(a1, a2) * radian_to_val) - theta_ai_max[1]);

	a1 = YbiOxy(2, y, phi) - ya[2];
	a2 = XbiOxy(2, x, phi) - xa[2];
	func_values.push_back((atan2(a1, a2) * radian_to_val) - theta_ai_max[2]);
}

//------------------------------------------------------------------------------------------
void g2i_6(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double a1 = YbiOxy(0, y, phi) - ya[0];
	double a2 = XbiOxy(0, x, phi) - xa[0];
	func_values.push_back(theta_ai_min[0] - (atan2(a1, a2) * radian_to_val));

	a1 = YbiOxy(1, y, phi) - ya[1];
	a2 = XbiOxy(1, x, phi) - xa[1];
	func_values.push_back(theta_ai_min[1] - (atan2(a1, a2) * radian_to_val));

	a1 = YbiOxy(2, y, phi) - ya[2];
	a2 = XbiOxy(2, x, phi) - xa[2];
	func_values.push_back(theta_ai_min[2] - (atan2(a1, a2) * radian_to_val));
}

//------------------------------------------------------------------------------------------
void g2i_11(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double diff_x = XbiOxy(0, x, phi) - xa[0];
	double diff_y = YbiOxy(0, y, phi) - ya[0];
	double rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	double theta_ai;

	if (diff_y >= 0.0)
		theta_ai = acos(diff_x / rho_i) * radian_to_val;
	else
		theta_ai = TWO_PI - (acos(diff_x / rho_i) * radian_to_val);

	if ((theta_ai - phi + M_PI) < TWO_PI)
		func_values.push_back( theta_ai - phi + M_PI - theta_bi_max[0]);
	else
		func_values.push_back(theta_ai - phi - M_PI - theta_bi_max[0]);


	diff_x = XbiOxy(1, x, phi) - xa[1];
	diff_y = YbiOxy(1, y, phi) - ya[1];
	rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	if (diff_y >= 0.0)
		theta_ai = acos(diff_x / rho_i) * radian_to_val;
	else
		theta_ai = TWO_PI - (acos(diff_x / rho_i) * radian_to_val);

	if ((theta_ai - phi + M_PI) < TWO_PI)
		func_values.push_back(theta_ai - phi + M_PI - theta_bi_max[1]);
	else
		func_values.push_back(theta_ai - phi - M_PI - theta_bi_max[1]);


	diff_x = XbiOxy(2, x, phi) - xa[2];
	diff_y = YbiOxy(2, y, phi) - ya[2];
	rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	if (diff_y >= 0.0)
		theta_ai = acos(diff_x / rho_i) * radian_to_val;
	else
		theta_ai = TWO_PI - (acos(diff_x / rho_i) * radian_to_val);

	if ((theta_ai - phi + M_PI) < TWO_PI)
		func_values.push_back(theta_ai - phi + M_PI - theta_bi_max[2]);
	else
		func_values.push_back(theta_ai - phi - M_PI - theta_bi_max[2]);
}

//------------------------------------------------------------------------------------------
void g2i_12(const double& x, const double& y, const double& phi, vector<double>& func_values)
{
	double diff_x = XbiOxy(0, x, phi) - xa[0];
	double diff_y = YbiOxy(0, y, phi) - ya[0];
	double rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	double theta_ai;

	if (diff_y >= 0.0)
		theta_ai = acos(diff_x / rho_i) * radian_to_val;
	else
		theta_ai = TWO_PI - (acos(diff_x / rho_i) * radian_to_val);

	if ((theta_ai - phi + M_PI) < TWO_PI)
		func_values.push_back(theta_bi_min[0] - theta_ai + phi - M_PI);
	else
		func_values.push_back(theta_bi_min[0] - theta_ai + phi + M_PI);


	diff_x = XbiOxy(1, x, phi) - xa[1];
	diff_y = YbiOxy(1, y, phi) - ya[1];
	rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	if (diff_y >= 0.0)
		theta_ai = acos(diff_x / rho_i) * radian_to_val;
	else
		theta_ai = TWO_PI - (acos(diff_x / rho_i) * radian_to_val);

	if ((theta_ai - phi + M_PI) < TWO_PI)
		func_values.push_back(theta_bi_min[1] - theta_ai + phi - M_PI);
	else
		func_values.push_back(theta_bi_min[1] - theta_ai + phi + M_PI);

	diff_x = XbiOxy(2, x, phi) - xa[2];
	diff_y = YbiOxy(2, y, phi) - ya[2];
	rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	if (diff_y >= 0.0)
		theta_ai = acos(diff_x / rho_i) * radian_to_val;
	else
		theta_ai = TWO_PI - (acos(diff_x / rho_i) * radian_to_val);

	if ((theta_ai - phi + M_PI) < TWO_PI)
		func_values.push_back(theta_bi_min[2] - theta_ai + phi - M_PI);
	else
		func_values.push_back(theta_bi_min[2] - theta_ai + phi + M_PI);
}

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
		if (func_values.back() < 0.0)				// func_values.back() - максимум
			solution.push_back(box);				// входит во множество решений
		else if (func_values.front() > 0.0)			// func_values.front() - минимум
			not_solution.push_back(box);			// в таком случае рассматриваемый box не входит во множество решений
		else										// подлежит разбинению и дальнейшему анализу
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
	func_values.reserve(points_in_3d);

	GetFusncValues(box, func_values);

	//printf("max_values vector size = %Iu\n", max_values.size());

	ClasifyBox(box, func_values);
}

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
	// вычитаем единицу т.к. кол-во отрезков, на которые делится ось координат 
	// всегда на единицу меньше, чем кол-во точек на этой оси
	double step_y = y_range / parts_per_axis;

	// определение границ по оси phi
	double phi_min = box_params[4];
	double phi_range = box_params[5];
	double phi_max = phi_min + phi_range;
	// определение шага сетки по координате phi
	double step_phi = phi_range / parts_per_axis;

	// координаты первой рассматриваемой точки на сетке
	double x = x_min, y = y_min, phi = phi_min;
	
	// находим значение всех 18 функций в начальной точке
	//GetFuncValuesInPoint(x, y, phi, gj_values);
	
	//cilk_for(unsigned k = 1; k <= total_amount_of_partions; ++k)

	// цикл по всем точкам, лежащим внутри box-а
	for (unsigned k = 0; k < points_in_3d; ++k)
	{ 
		// определяем координаты очередной точки на сетке
		x = x_min + step_x*(k % points_per_axis);
		y = y_min + step_y*((k / points_per_axis) % points_per_axis);
		phi = phi_min + step_phi*(k / squared_points_per_axis);
		
		GetFuncValuesInPoint(x, y, phi, gj_values);
	}
	
	std::sort(gj_values.begin(), gj_values.end());
}

high_level_analysis::high_level_analysis(const double* init_params) :
	low_level_fragmentation(init_params) {}

high_level_analysis::high_level_analysis(Box& box) : low_level_fragmentation(box) {}

/// 
void high_level_analysis::GetFusncValues(const Box& box, vector<double>& max_vals)
{
	//puts("GetMaxs()");

	double curr_box_diagonal = box.GetDiagonal();
	//printf("cur_box diagonal = %lf\n", curr_box_diagonal);

	if (curr_box_diagonal > g_precision)
	{
		double* box_params = new double[6];
		box.GetParameters(box_params);

		// нахождение максимумов всех 18 функций в рассматриваемом box-е
		GetGjValues(box_params, max_vals);
	}
	else	// достигнута желаемая точность аппроксимации рабочего пространства
	{
		boundary.push_back(box);
	}
}

///
void high_level_analysis::GetSolution()
{
	// initial box is defined in main()
//	vector<double> initial_box_params = { -20.0, 40.0, -20.0, 40.0, 0.0, 360.0 };
//	current_box = Box(initial_box_params);

	unsigned int tree_depth = FindTreeDepth();		// find out number of outer loop iterations

	printf("Tree depth is %d \n", tree_depth);

	//next_iter_boxes->push_back(current_box);
	next_iter_boxes.push_back(current_box);

	vector<Box> curr_iter_boxes, empty_vec;

	puts("Getting solution");

	for (unsigned int level = 0; level < (tree_depth + 1); ++level)
	{
		//printf("Level = %d\n", level);

		size_t number_of_box_on_level = next_iter_boxes.size();

		//next_iter_boxes.move_out(curr_iter_boxes);
		curr_iter_boxes = next_iter_boxes;

		// ! необходимо очищать содержимое next_iter_boxes
		//next_iter_boxes.set_value(empty_vec);
		next_iter_boxes.clear();

		//printf("level is %d ",level);

		//cilk_for(size_t i = 0; i < number_of_box_on_level; ++i)
		for (size_t i = 0; i < number_of_box_on_level; ++i)
		{
			//printf("i is %d\n", i);
			GetBoxType(curr_iter_boxes[i]);
		}
	}
}






///// 
//void WriteResults( const char* file_names[], unsigned int number_of_files )
//{
//	if (number_of_files != 3)
//	{
//		puts("NOT 3 elments!");
//		return;
//	}
//
//	std::ofstream fout(file_names[0]);
//
//	if (fout.is_open() == false)
//	{
//		puts("Can't open file for writing!");
//		return;
//	}
//	else
//	{
//		//vector<Box> result_solution;
//		//solution.move_out(result_solution);
//
//		//double left, bottom, width, height;
//
//		std::vector<Box>::const_iterator res_it_beg = solution.begin();
//		std::vector<Box>::const_iterator res_it_end = solution.end();
//
//		vector<double> box_params;
//		box_params.reserve(6);
//
//		for (res_it_beg; res_it_beg != res_it_end; ++res_it_beg)
//		{
//			res_it_beg->GetParameters(box_params);
//
//			fout << box_params[0] << ' ' << box_params[2] << ' ' << box_params[1] << ' ' << box_params[3] << '\n';
//		}
//
//		fout.close();
//	}
//
//
//	std::ofstream fout2(file_names[1]);
//
//	if (fout2.is_open() == false)
//	{
//		puts("Ошибка при открытии файла!");
//		return;
//	}
//	else
//	{
//		//vector<Box> result_solution;
//		//solution.move_out(result_solution);
//
//		//double left, bottom, width, height;
//
//		std::vector<Box>::const_iterator res_it_beg = not_solution.begin();
//		std::vector<Box>::const_iterator res_it_end = not_solution.end();
//
//		vector<double> box_params;
//		box_params.reserve(6);
//
//		for (res_it_beg; res_it_beg != res_it_end; ++res_it_beg)
//		{
//			res_it_beg->GetParameters(box_params);
//
//			fout2 << box_params[0] << ' ' << box_params[2] << ' ' << box_params[1] << ' ' << box_params[3] << '\n';
//		}
//
//		fout2.close();
//	}
//
//	std::ofstream fout3(file_names[2]);
//
//	if (fout3.is_open() == false)
//	{
//		puts("Ошибка при открытии файла!");
//		return;
//	}
//	else
//	{
//		//vector<Box> result_solution;
//		//solution.move_out(result_solution);
//
//		//double left, bottom, width, height;
//
//		std::vector<Box>::const_iterator res_it_beg = not_solution.begin();
//		std::vector<Box>::const_iterator res_it_end = not_solution.end();
//
//		vector<double> box_params;
//		box_params.reserve(6);
//
//		for (res_it_beg; res_it_beg != res_it_end; ++res_it_beg)
//		{
//			res_it_beg->GetParameters(box_params);
//
//			fout3 << box_params[0] << ' ' << box_params[2] << ' ' << box_params[1] << ' ' << box_params[3] << '\n';
//		}
//
//		fout3.close();
//	}
//}


void PrintWorkspace()
{
	puts("Printing workspace!");

	Engine* matlab_ptr;
	matlab_ptr = engOpen("null");

	mxArray* matlab_x_min = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxArray* matlab_x_max = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxArray* matlab_y_min = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxArray* matlab_y_max = mxCreateDoubleMatrix(1, 1, mxREAL);

	double* vs_x_min = mxGetPr(matlab_x_min);
	double* vs_x_max = mxGetPr(matlab_x_max);
	double* vs_y_min = mxGetPr(matlab_y_min);
	double* vs_y_max = mxGetPr(matlab_y_max);

	engEvalString(matlab_ptr, "figure('units','normalized','outerposition',[0 0 1 1]),");
	engEvalString(matlab_ptr, "hold on,");

	size_t sol_sz = solution.size();
	size_t not_sol_sz = not_solution.size();
	size_t boundary_sz = boundary.size();

	printf("Size of solution vector = %Iu\n", sol_sz);
		
	printf("Size of not_solution vector = %Iu\n", not_sol_sz);
	printf("Size of boundary vector = %Iu\n", boundary_sz);

	std::vector<Box>::const_iterator solution_it_begin = solution.begin();
	std::vector<Box>::const_iterator solution_it_end = solution.end();

	double* box_params = new double[6];
	double phi_min, phi_max;

	for (solution_it_begin; solution_it_begin != solution_it_end; ++solution_it_begin)
	{

		solution_it_begin->GetParameters(box_params);

		*vs_x_min = box_params[0];
		*vs_x_max = *vs_x_min + box_params[1];

		*vs_y_min = box_params[2];
		*vs_y_max = *vs_y_min + box_params[3];

		phi_min = box_params[4];
		phi_max = phi_min + box_params[5];

		printf("phi_min = %lf  phi_max = %lf\n", phi_min, phi_max);

		if (cmp_angle >= phi_min && cmp_angle <= phi_max)
		{
			puts("plotting!");

			engPutVariable(matlab_ptr, "x_min", matlab_x_min);
			engPutVariable(matlab_ptr, "x_max", matlab_x_max);
			engPutVariable(matlab_ptr, "y_min", matlab_y_min);
			engPutVariable(matlab_ptr, "y_max", matlab_y_max);

			engEvalString(matlab_ptr, "line([x_min x_max], [y_min, y_min], 'color', 'green'),");
			engEvalString(matlab_ptr, "line([x_min x_min], [y_min, y_max], 'color', 'green'),");
			engEvalString(matlab_ptr, "line([x_min x_max], [y_max, y_max], 'color', 'green'),");
			engEvalString(matlab_ptr, "line([x_max x_max], [y_min, y_max], 'color', 'green')");
		}
	}

	////engEvalString(matlab_ptr, "axis([0 8 0 8])");

	//system("pause");

	// free memory
	/*mxDestroyArray(matlab_x_min);
	mxDestroyArray(matlab_x_max);
	mxDestroyArray(matlab_y_min);
	mxDestroyArray(matlab_y_max);

	engClose(matlab_ptr);*/
}