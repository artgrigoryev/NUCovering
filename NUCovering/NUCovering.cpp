#define _USE_MATH_DEFINES
#include <cmath>
#include "fragmentation.h"

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
const double g_precision = 0.25;

/// количество отрезков, на которые разделяется грань box-а 
const unsigned int number_of_partitions_per_axis = 64;


int main()
{
	//	setlocale(LC_ALL,"Rus");

	//	__cilkrts_end_cilk();
	high_level_analysis main_object;
	//	__cilkrts_set_param("nworkers", "1");
	//	cout << "Number of workers " << __cilkrts_get_nworkers() << endl;
	//	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	main_object.GetSolution();
	//	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	//	cout << "Duration is: " << time_span.count() << " seconds" << endl;
	//	auto duration = duration_cast<microseconds>(t2 - t1).count();
	//	cout << duration << " sec";

	return 0;
}