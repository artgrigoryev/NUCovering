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

/// ���������� ����� �� ������ �� ����
extern const unsigned int points_per_axis;

//extern const double cmp_angle;

typedef std::pair<Box, Box> boxes_pair;


class low_level_fragmentation
{
protected:
	Box current_box;
	unsigned int FindTreeDepth();
	void SplitByX(const Box& box, boxes_pair& new_pair_of_boxes) const;		// ��������� box'� �� ���������� �
	void SplitByY(const Box& box, boxes_pair& new_pair_of_boxes) const;		// ��������� box'� �� ���������� y
	void SplitByPhi(const Box& box, boxes_pair& new_pair_of_boxes) const;	// ��������� box'a �� ���������� z
	void GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes) const;				// ��������� ���� ����� box'�� ����� ���������
	void ClasifyBox(const Box& box, const std::vector<double>& func_values) const;				// ������� ������� box'�, ����������� ���� �������������� box'�
	void GetBoxType(const Box& box);
	virtual void GetFusncValues(const Box& box, std::vector<double>& max_vals) = 0;		// ��������� ��������� � ���������� ��� ������� gj

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
