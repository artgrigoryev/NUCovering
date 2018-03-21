//#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <chrono>
#include "fragmentation.h"

using namespace std::chrono;
using std::vector;

/// ���������� ����� �[1,2,3] � ������� ��������� Oxy
const double xa[] = { -15.0, 15.0, 0.0 };
const double ya[] = { -5.0 * sqrt(3.0), -5.0 * sqrt(3.0), 10.0 * sqrt(3.0) };

/// ���������� ����� B[1,2,3] � ������� ��������� O'xy
const double xb[] = { -5.0, 5.0, 0.0 };
const double yb[] = { -5.0 * sqrt(3.0) / 3.0, -5.0 * sqrt(3.0) / 3.0, 10.0 * sqrt(3.0) / 3.0 };

/// ����������� ����� �����, �� ������� ��������� ����������� ��������
const double p_min = 12.0;
/// ������������ ����� �����, �� ������� ��������� ����������� ��������
const double p_max = 27.0;

/// ����������� ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� A[1,2,3] ������������ ������� ��������� Oxy
const double theta_ai_min[] = { 10.0, 0.0, 0.0 };
/// ������������ ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� A[1,2,3] ������������ ������� ��������� Oxy
const double theta_ai_max[] = { 350.0, 300.0, 360.0 };

/// ����������� ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� B[1,2,3] ������������ ����������� ������� ��������� O'xy
const double theta_bi_min[] = { 10.0, 30.0, 0.0 };
/// ������������ ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� B[1,2,3] ������������ ����������� ������� ��������� O'xy
const double theta_bi_max[] = { 340.0, 360.0, 330.0 };


/// �������� ������������� �������� ������������
const double g_precision = 0.8;

/// ���������� ����� �� ������ �� ���� 
const unsigned int points_per_axis = 12;


#pragma comment (lib, "libmat.lib")
#pragma comment (lib, "libmx.lib")
#pragma comment (lib, "libmex.lib")
#pragma comment (lib, "libeng.lib")


int main()
{
	//	__cilkrts_end_cilk();

	// initial box parameters
	double initial_box_params[] = { -20.0, 40.0, -20.0, 40.0, 0.0, 360.0 };
	high_level_analysis main_object(initial_box_params);

	//	__cilkrts_set_param("nworkers", "1");
	//	cout << "Number of workers " << __cilkrts_get_nworkers() << endl;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	main_object.GetSolution();

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	auto duration = duration_cast<seconds>(t2 - t1).count();
	std::cout << std::endl << "Duration is: " << duration << " seconds" << std::endl;
	std::cout << "Duration is: " << duration/60 << " minutes " << duration%60 << "seconds" << std::endl << std::endl;

	// ������ �����, ������������ ������� ������ �������� ����������� �������� ������������
	double cmp_angles[] = { 50.0, 80.0, 120.0, 140.0 };

	PrintWorkspace(cmp_angles, sizeof(cmp_angles) / sizeof(cmp_angles[0]));

	return 0;
}