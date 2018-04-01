#pragma once
#include "box.h"

/// координаты точек А[1,2,3] в системе координат Oxy
extern const double xa[];
extern const double ya[];

/// координаты точек B[1,2,3] в системе координат O'xy
extern const double xb[];
extern const double yb[];

/// минимальная длина штанг, на которые крепиться треугольная площадка
extern const double p_min;
/// максимальная длина штанг, на которые крепиться треугольная площадка
extern const double p_max;

/// минимальный угол, на который может повернутся штанга 
/// в точке закрепления A[1,2,3] относительно системы координат Oxy
extern const double theta_ai_min[];
/// максимальный угол, на который может повернутся штанга 
/// в точке закрепления A[1,2,3] относительно системы координат Oxy
extern const double theta_ai_max[];

/// минимальный угол, на который может повернутся штанга 
/// в точке закрепления B[1,2,3] относительно вращающейся системы координат O'xy
extern const double theta_bi_min[];
/// максимальный угол, на который может повернутся штанга 
/// в точке закрепления B[1,2,3] относительно вращающейся системы координат O'xy
extern const double theta_bi_max[];

/// точность аппроксимации рабочего пространства
extern const double g_precision;

/// количество точек на каждой из осей
extern const unsigned int points_per_axis;

//extern const double cmp_angle;

typedef std::pair<Box, Box> boxes_pair;


class low_level_fragmentation
{
protected:
	Box current_box;
	unsigned int FindTreeDepth();
	void SplitByX(const Box& box, boxes_pair& new_pair_of_boxes) const;		// разбиение box'а по координате х
	void SplitByY(const Box& box, boxes_pair& new_pair_of_boxes) const;		// разбиение box'а по координате y
	void SplitByPhi(const Box& box, boxes_pair& new_pair_of_boxes) const;	// разбиение box'a по координате z
	void GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes) const;				// получение двух новых box'ов путем разбиения
	void ClasifyBox(const Box& box, const std::vector<double>& func_values) const;				// функция анализа box'а, определения типа принадлежности box'а
	void GetBoxType(const Box& box);
	virtual void GetFusncValues(const Box& box, std::vector<double>& max_vals) = 0;		// получение минимумов и максимумов для функций gj

public:
	/// default constructor
	low_level_fragmentation() {}

	/// parametrized constructor
	low_level_fragmentation(const double* init_params);

	/// copy constructor
	low_level_fragmentation(const Box &box);

	/// virtual destructor
	virtual ~low_level_fragmentation() {}
};


class high_level_analysis : public low_level_fragmentation
{
protected:
	void GetFusncValues(const Box& box, std::vector<double>& max_vals) override;
public:
	/// default constructor
	high_level_analysis() {}

	/// parametrized constructor
	high_level_analysis(const double* init_params);

	/// copy constructor
	high_level_analysis(Box& box);

	void GetSolution();
};

void PrintWorkspace(const double* cmp_angles, const unsigned& cmp_angles_size,
                     const unsigned& rings_count, double& phi_min, double& phi_max);
