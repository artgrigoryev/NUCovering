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

/// количество отрезков, на которые разделяется грань box-а 
extern const unsigned int number_of_partitions_per_axis;


typedef std::pair<Box, Box> boxes_pair;


class low_level_fragmentation
{
protected:
	Box current_box;
	unsigned int FindTreeDepth();
	void VerticalSplitter(const Box& box, boxes_pair& vertical_splitter_pair) const;		// разбиение box'а по длине
	void HorizontalSplitter(const Box& box, boxes_pair& horizontal_splitter_pair) const;	// разбиение box'а по высоте
	void GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes) const;				// получение двух новых box'ов путем разбиения
	int ClasifyBox(const std::vector<double>& maxs_vec) const;				// функция анализа box'а, определения типа принадлежности box'а
	void GetBoxType(const Box& box);
	virtual void GetMaxs(const Box& box, std::vector<double>& max_vals) = 0;		// получение минимумов и максимумов для функций gj

public:
	/// default constructor
	low_level_fragmentation() {}

	/// parametrized constructor
	low_level_fragmentation(const std::vector<double>& init_vec);

	/// copy constructor
	low_level_fragmentation(const Box &box);

	/// virtual destructor
	virtual ~low_level_fragmentation() {}
};


class high_level_analysis : public low_level_fragmentation
{
protected:
	void GetMaxs(const Box& box, std::vector<double>& max_vals);
public:
	/// default constructor
	high_level_analysis() {}

	/// parametrized constructor
	high_level_analysis(std::vector<double>& init_vec);

	/// copy constructor
	high_level_analysis(Box& box);

	void GetSolution();

	//void WriteResults( const char* file_names[] );
};