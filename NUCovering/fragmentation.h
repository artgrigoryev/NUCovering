#pragma once
#include "box.h"

/// ���������� ����� �[1,2,3] � ������� ��������� Oxy
extern const double xa[];
extern const double ya[];

/// ���������� ����� B[1,2,3] � ������� ��������� O'xy
extern const double xb[];
extern const double yb[];

/// ����������� ����� �����, �� ������� ��������� ����������� ��������
extern const double p_min;
/// ������������ ����� �����, �� ������� ��������� ����������� ��������
extern const double p_max;

/// ����������� ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� A[1,2,3] ������������ ������� ��������� Oxy
extern const double theta_ai_min[];
/// ������������ ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� A[1,2,3] ������������ ������� ��������� Oxy
extern const double theta_ai_max[];

/// ����������� ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� B[1,2,3] ������������ ����������� ������� ��������� O'xy
extern const double theta_bi_min[];
/// ������������ ����, �� ������� ����� ���������� ������ 
/// � ����� ����������� B[1,2,3] ������������ ����������� ������� ��������� O'xy
extern const double theta_bi_max[];

/// �������� ������������� �������� ������������
extern const double g_precision;

/// ���������� ��������, �� ������� ����������� ����� box-� 
extern const unsigned int number_of_partitions_per_axis;


typedef std::pair<Box, Box> boxes_pair;


class low_level_fragmentation
{
protected:
	Box current_box;
	unsigned int FindTreeDepth();
	void VerticalSplitter(const Box& box, boxes_pair& vertical_splitter_pair) const;		// ��������� box'� �� �����
	void HorizontalSplitter(const Box& box, boxes_pair& horizontal_splitter_pair) const;	// ��������� box'� �� ������
	void GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes) const;				// ��������� ���� ����� box'�� ����� ���������
	int ClasifyBox(const std::vector<double>& maxs_vec) const;				// ������� ������� box'�, ����������� ���� �������������� box'�
	void GetBoxType(const Box& box);
	virtual void GetMaxs(const Box& box, std::vector<double>& max_vals) = 0;		// ��������� ��������� � ���������� ��� ������� gj

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