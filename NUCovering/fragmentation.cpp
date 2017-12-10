#define _USE_MATH_DEFINES
#include <cmath>
#include <cilk/cilk.h>
#include <cilk/reducer_vector.h>
#include "fragmentation.h"

using std::vector;

const double p_min2 = p_min*p_min;
const double p_max2 = p_max*p_max;
const double two_pi = M_PI * 2;

// общее количество отрезков для 3-х осей координат
const unsigned int total_amount_of_partions = pow(number_of_partitions_per_axis, 3);

const unsigned int partions2 = number_of_partitions_per_axis*number_of_partitions_per_axis;

// количество функций gj() в группе gj[]
const unsigned int funcs_per_group = 3;
// коэффициент, используемый для перевода значения угла из градусов в радианы
const double koeff = M_PI / 180;

//------------------------------------------------------------------------------------------
double Xci(const unsigned int& i, const double& phi)
{
	return xa[i] - xb[i] * cos(phi * koeff) + yb[i] * sin(phi * koeff);
}

//------------------------------------------------------------------------------------------
double Yci(const unsigned int& i, const double& phi)
{
	return ya[i] - xb[i] * sin(phi * koeff) - yb[i] * cos(phi* koeff);
}

//------------------------------------------------------------------------------------------
// находим координату x произвольной точки, лежащей во
// вращающейся плоскости O'xy, относительно неподвижной плоскости Oxy
double XbiOxy(const unsigned int& i, const double& x, const double& phi)
{
	return x + xb[i] * cos(phi * koeff) - yb[i] * sin(phi * koeff);
}

//------------------------------------------------------------------------------------------
// находим координату y произвольной точки, лежащей во
// вращающейся плоскости O'xy, относительно неподвижной плоскости Oxy
double YbiOxy(const unsigned int& i, const double& y, const double& phi)
{
	return y + xb[i] * sin(phi* koeff) + yb[i] * cos(phi * koeff);
}

/// functions gj()
//------------------------------------------------------------------------------------------
double g2i_1(const unsigned int& i, const double& x, const double& y, const double& phi)
{
	// params[0] -> x, param[1] -> y, params[2] -> phi
	double a1 = XbiOxy(i, x, y) - xa[i];
	double a2 = YbiOxy(i, y, phi) - ya[i];
	return a1*a1 + a2*a2 - p_max2;
}

//------------------------------------------------------------------------------------------
double g2i(const unsigned int& i, const double& x, const double& y, const double& phi)
{
	// params[0] -> x, param[1] -> y, params[2] -> phi
	double a1 = XbiOxy(i, x, phi) - Xci(i, phi);
	double a2 = YbiOxy(i, y, phi) - Yci(i, phi);
	return p_min2 - a1*a1 - a2*a2;
}

//------------------------------------------------------------------------------------------
double g2i_5(const unsigned int& i, const double& x, const double& y, const double& phi)
{
	// params[0] -> x, param[1] -> y, params[2] -> phi
	double a1 = YbiOxy(i, y, phi) - ya[i];
	double a2 = XbiOxy(i, x, phi) - xa[i];
	return atan2(a1, a2) - theta_ai_max[i];
}

//------------------------------------------------------------------------------------------
double g2i_6(const unsigned int& i, const double& x, const double& y, const double& phi)
{
	// params[0] -> x, param[1] -> y, params[2] -> phi
	double a1 = YbiOxy(i, y, phi) - ya[i];
	double a2 = XbiOxy(i, x, phi) - xa[i];
	return theta_ai_min[i] - atan2(a1, a2);
}

//------------------------------------------------------------------------------------------
double g2i_11(const unsigned int& i, const double& x, const double& y, const double& phi)
{
	// params[0] -> x, param[1] -> y, params[2] -> phi
	double diff_x = XbiOxy(i, x, phi) - xa[i];
	double diff_y = YbiOxy(i, y, phi) - ya[i];
	double rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	double theta_ai;

	if (diff_y >= 0)
		theta_ai = acos(diff_x / rho_i);
	else
		theta_ai = two_pi - acos(diff_x / rho_i);

	if ((theta_ai - phi + M_PI) < two_pi)
		return theta_ai - phi + M_PI - theta_bi_max[i];
	else
		return theta_ai - phi - M_PI - theta_bi_max[i];
}

//------------------------------------------------------------------------------------------
double g2i_12(const unsigned int& i, const double& x, const double& y, const double& phi)
{
	// params[0]->x, param[1]->y, params[2]->phi
	double diff_x = XbiOxy(i, x, phi) - xa[i];
	double diff_y = YbiOxy(i, y, phi) - ya[i];
	double rho_i = sqrt(diff_x * diff_x + diff_y * diff_y);

	double theta_ai;

	if (diff_y >= 0)
		theta_ai = acos(diff_x / rho_i);
	else
		theta_ai = two_pi - acos(diff_x / rho_i);

	if ((theta_ai - phi + M_PI) < two_pi)
		return theta_bi_min[i] - theta_ai + phi - M_PI;
	else
		return theta_bi_min[i] - theta_ai + phi + M_PI;
}


// f_ptr - указатель на группу функций gj()
typedef double(*f_ptr) (const unsigned int& i, const double& x, const double& y, const double& phi);

// массив указателей на группы функций gj()
f_ptr gj[] = { g2i_1, g2i, g2i_5, g2i_6, g2i_11, g2i_12 };


/// контейнер box-ов, входящих в аппроксимацию рабочего простраства
cilk::reducer< cilk::op_vector<Box> > solution;
/// контейнер box-ов, не входящих в аппроксимацию рабочего простраства
cilk::reducer< cilk::op_vector<Box> > not_solution;
/// контейнер граничных box-ов, которые не относятся ни к solution, ни к not_solution
cilk::reducer< cilk::op_vector<Box> > boundary;
/// контейнер box-ов, обрабатываемых на следующей итерации алгоритма
cilk::reducer< cilk::op_vector<Box> > next_iter_boxes;


//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(const vector<double>& init_vec)
{
	current_box = Box(init_vec);
}

//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(const Box& box)
{
	current_box = box;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::VerticalSplitter(const Box& box, boxes_pair& vertical_splitter_pair) const
{
	vector<double> box_params;
	box.GetParameters(box_params);

	double x_min = box_params[0];
	double new_widht = box_params[1] * 0.5;
	double new_angle_min = box_params[4];
	double new_angle_range = box_params[5] * 0.5;

	vector<double> new_box_params = { x_min, new_widht, box_params[2], box_params[3], new_angle_min, new_angle_range };
	vertical_splitter_pair.first = Box(new_box_params);

	new_box_params[0] = x_min + new_widht;
	new_box_params[4] = new_angle_min + new_angle_range;
	vertical_splitter_pair.second = Box(new_box_params);
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::HorizontalSplitter(const Box& box, boxes_pair& horizontal_splitter_pair) const
{
	vector<double> box_params;
	box.GetParameters(box_params);

	double y_min = box_params[2];
	double new_height = box_params[3] * 0.5;
	double new_angle_min = box_params[4];
	double new_angle_range = box_params[5] * 0.5;

	vector<double> new_box_params = { box_params[0], box_params[1], y_min, new_height, new_angle_min, new_angle_range };
	horizontal_splitter_pair.first = Box(new_box_params);

	new_box_params[2] = y_min + new_height;
	new_box_params[4] = new_angle_min + new_angle_range;
	horizontal_splitter_pair.second = Box(new_box_params);
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes) const
{
	double param_width, param_height;
	box.GetWidhtHeight(param_width, param_height);

	if (param_width >= param_height)
		VerticalSplitter(box, new_pair_of_boxes);
	else
		HorizontalSplitter(box, new_pair_of_boxes);
}

//------------------------------------------------------------------------------------------
unsigned int low_level_fragmentation::FindTreeDepth()
{
	double box_diagonal = current_box.GetDiagonal();

	if (box_diagonal <= g_precision)
	{
		return 0;
	}
	else
	{
		boxes_pair new_boxes;
		VerticalSplitter(current_box, new_boxes);
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
int low_level_fragmentation::ClasifyBox(const std::vector<double>& max_values) const
{
	std::vector<double>::const_iterator maxs_iter_beg = max_values.begin();
	std::vector<double>::const_iterator maxs_iter_end = max_values.end();

	double min_max = *maxs_iter_beg;
	double max_max = *maxs_iter_beg;

	for (maxs_iter_beg; maxs_iter_beg != maxs_iter_end; ++maxs_iter_beg)
	{
		if (*maxs_iter_beg < min_max)			// проверка условия min gj(x) > 0 для любого j из 1,...,m
		{
			min_max = *maxs_iter_beg;
		}
		if (*maxs_iter_beg > max_max)
		{
			max_max = *maxs_iter_beg;
		}
	}

	if (max_max < 0)
		return 0;							// входит во множество решений
	else if (min_max > 0)
		return 1;							// в таком случае рассматриваемый box не входит во множество решений
	else
		return 2;							// подлежит разбинению и дальнейшему анализу
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetBoxType(const Box& box)
{
	// вектор хранящий значения максимумов функций gj[] в пределах box'a 
	vector<double> max_values;

	GetMaxs(box, max_values);

	int res = ClasifyBox(max_values);

	if (res == 0)
	{
		solution->push_back(box);					// входит во множество решений
	}
	else if (res == 1)
	{
		not_solution->push_back(box);					// в таком случае рассматриваемый box не входит во множество решений
	}
	else
	{
		boxes_pair another_boxes;
		GetNewBoxes(box, another_boxes);				// подлежит разбинению и дальнейшему анализу

		next_iter_boxes->push_back(another_boxes.first);
		next_iter_boxes->push_back(another_boxes.second);
	}
}

void GetFuncsMaxValues(const vector<double>& box_params, vector<double>& gj_max_vals,
	double(*func)(const unsigned int&, const double&, const double&, const double&))
{
	/// определение координат точки, принадлежащей box-у
	// границы по оси x
	double x_min = box_params[0];
	double x_max = x_min + box_params[1];
	// границы по оси y
	double y_min = box_params[2];
	double y_max = y_min + box_params[3];
	// границы по оси phi
	double phi_min = box_params[4];
	double phi_max = phi_min + box_params[5];

	// определение шага по x
	double step_x = (x_min - x_max) / number_of_partitions_per_axis;
	// определение шага по y
	double step_y = (y_max - y_min) / number_of_partitions_per_axis;
	// определение шага по phi
	double step_phi = (phi_max - phi_min) / number_of_partitions_per_axis;

	// координаты первой рассматриваемой точки
	double x = x_min, y = y_min, phi = phi_min;
	double max_funcs_values[funcs_per_group];

	for (unsigned int i = 0; i < funcs_per_group; ++i)
	{
		// считается, что максимумы функций gj[] достигаются в первой точке
		max_funcs_values[i] = func(i, x, y, phi);
	}

	double funcs_values[funcs_per_group];

	// цикл по всем точкам, лежащим внутри box-а
	cilk_for(unsigned int k = 1; k <= total_amount_of_partions; ++k)
	{
		x = x_min + step_x*(k % number_of_partitions_per_axis);
		y = y_min + step_y*((k / number_of_partitions_per_axis) % number_of_partitions_per_axis);
		phi = phi_min + k*((k / partions2) % number_of_partitions_per_axis);

		// определение значения 3-х функций gj[] в новой точке и
		// нахождение "текущего" максимума этих 3-х функций
		for (unsigned int i = 0; i < funcs_per_group; ++i)
		{
			funcs_values[i] = func(i, x, y, phi);
			if (funcs_values[i] > max_funcs_values[i])
				max_funcs_values[i] = funcs_values[i];
		}
	}

	gj_max_vals.push_back(max_funcs_values[0]);
	gj_max_vals.push_back(max_funcs_values[1]);
	gj_max_vals.push_back(max_funcs_values[2]);
}

high_level_analysis::high_level_analysis(std::vector<double>& init_vec) :
	low_level_fragmentation(init_vec) {}

high_level_analysis::high_level_analysis(Box& box) : low_level_fragmentation(box) {}

/// 
void high_level_analysis::GetMaxs(const Box& box, vector<double>& max_vals)
{
	double curr_box_diagonal = box.GetDiagonal();

	if (curr_box_diagonal > g_precision)
	{
		vector<double> box_params;
		box_params.reserve(6);
		box.GetParameters(box_params);

		// нахождение максимумов всех 18 функций в рассматриваемом box-е
		for (size_t i = 0; i < sizeof(gj) / sizeof(gj[0]); ++i)
		{
			GetFuncsMaxValues(box_params, max_vals, gj[i]);
		}
	}
	else	// достигнута желаемая точность аппроксимации рабочего пространства
	{
		boundary->push_back(box);
	}
}

///
void high_level_analysis::GetSolution()
{
	//  TODO: определить initial box
	//current_box = Box( (-1)*g_l1_max, 0.0, g_l1_max + g_l0 + g_l2_max, __min( g_l1_max, g_l2_max ), phi );

	unsigned int tree_depth = FindTreeDepth();		// find out number of outer loop iterati

													//printf("Tree depth is %d \n", tree_depth);

	next_iter_boxes->push_back(current_box);

	vector<Box> curr_iter_boxes, empty_vec;

	for (unsigned int level = 0; level < (tree_depth + 1); ++level)
	{
		size_t number_of_box_on_level = (1 << level);

		next_iter_boxes.move_out(curr_iter_boxes);

		// ! необходимо очищать содержимое next_iter_boxes
		next_iter_boxes.set_value(empty_vec);

		//printf("level is %d ",level);

		cilk_for(size_t i = 0; i < number_of_box_on_level; ++i)
		{
			//printf("i is %d\n", i);
			GetBoxType(curr_iter_boxes[i]);
		}
	}
}

/*
/// TODO: необходимо дописать
void WriteResults( const char* file_names[] )
{
if ( sizeof( file_names ) / sizeof( file_names[0] ) != 3 )
{
puts( "Количество элементов массива file_names[] должно быть 3!" );
return;
}

std::ofstream fout( file_names[0] );

if ( fout.is_open() == false )
{
puts( "Ошибка при открытии файла!" );
return;
}
else
{
vector<Box> result_solution;
solution.move_out( result_solution );

double left, bottom, width, height, angle;

std::vector<Box>::const_iterator res_it_beg = result_solution.begin();
std::vector<Box>::const_iterator res_it_end = result_solution.end();

vector<double> box_params;
box_params.reserve( 6 );

for ( res_it_beg; res_it_beg != res_it_end; ++res_it_beg )
{
res_it_beg->GetParameters( box_params );

fout << left << ' ' << bottom << ' ' << width << ' ' << height << '\n';
}

fout.close();
}
}
*/