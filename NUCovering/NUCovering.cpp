#define _USE_MATH_DEFINES
#include <cmath>
#include "fragmentation.h"

using std::vector;

/// координаты точек А[1,2,3] в системе координат Oxy
const double xa[] = { -15.0, 15.0, 0.0 };
const double ya[] = { -5.0 * sqrt(3.0), -5.0 * sqrt(3.0), 10.0 * sqrt(3.0) };

/// координаты точек B[1,2,3] в системе координат O'xy
const double xb[] = { -5.0, 5.0, 0.0 };
const double yb[] = { -5.0 * sqrt(3.0) / 3.0, -5.0 * sqrt(3.0) / 3.0, 10.0 * sqrt(3.0) / 3.0 };

/// минимальная длина штанг, на которые крепиться треугольная площадка
const double p_min = 12.0;
/// максимальная длина штанг, на которые крепиться треугольная площадка
const double p_max = 27.0;

/// минимальный угол, на который может повернутся штанга 
/// в точке закрепления A[1,2,3] относительно системы координат Oxy
const double theta_ai_min[] = { 10.0, 0.0, 0.0 };
/// максимальный угол, на который может повернутся штанга 
/// в точке закрепления A[1,2,3] относительно системы координат Oxy
const double theta_ai_max[] = { 350.0, 300.0, 360.0 };

/// минимальный угол, на который может повернутся штанга 
/// в точке закрепления B[1,2,3] относительно вращающейся системы координат O'xy
const double theta_bi_min[] = { 10.0, 30.0, 0.0 };
/// максимальный угол, на который может повернутся штанга 
/// в точке закрепления B[1,2,3] относительно вращающейся системы координат O'xy
const double theta_bi_max[] = { 340.0, 360.0, 330.0 };

/// точность аппроксимации рабочего пространства
const double g_precision = 2.0;

/// количество точек на каждой из осей 
const unsigned int points_per_axis = 4;

// значение угла, относительно которого строим проекцию
const double cmp_angle = 22.25;

//extern vector<Box> solution;
//extern vector<Box> not_solution;
//extern vector<Box> boundary;
//extern vector<Box> next_iter_boxes;


#pragma comment (lib, "libmat.lib")
#pragma comment (lib, "libmx.lib")
#pragma comment (lib, "libmex.lib")
#pragma comment (lib, "libeng.lib")


int main()
{
	//	setlocale(LC_ALL,"Rus");

	//	__cilkrts_end_cilk();
	// initial box parameters
	double initial_box_params[] = { -20.0, 40.0, -20.0, 40.0, 0.0, 360.0 };
	high_level_analysis main_object(initial_box_params);
	//	__cilkrts_set_param("nworkers", "1");
	//	cout << "Number of workers " << __cilkrts_get_nworkers() << endl;
	//	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	main_object.GetSolution();
	//	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	//	cout << "Duration is: " << time_span.count() << " seconds" << endl;
	//	auto duration = duration_cast<microseconds>(t2 - t1).count();
	//	cout << duration << " sec";
	
	/*const char* out_files[] = {"C:/Users/Artyom/Documents/Visual Studio 2015/Projects/NUCovering/solution.txt",
		"C:/Users/Artyom/Documents/Visual Studio 2015/Projects/NUCovering/not_solution.txt",
		"C:/Users/Artyom/Documents/Visual Studio 2015/Projects/NUCovering/boundary.txt" };

	WriteResults(out_files, 3);*/

	puts("Workspace founded!");

	PrintWorkspace();

	return 0;
}