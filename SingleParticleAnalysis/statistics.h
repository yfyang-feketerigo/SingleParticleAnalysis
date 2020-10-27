/*

class Statistics :ͳ���࣬����ͳ��һά���ݣ��ɼ��㣺
	���ֵ��Сֵ����ֵ����λ�����ֲ�

	struct Distribution���ֲ��ṹ������ֲ��������͡�
*/
#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>	
struct Distribution
{
	double min = 0; //�����е���Сֵ
	double max = 0; //�����е����ֵ
	std::vector<double> distribution; //�ֲ�
	void normalization() //��һ��
	{
		double sum = 0;
		for (size_t i = 0; i < distribution.size(); i++)
		{
			sum += distribution[i];
		}
		for (size_t i = 0; i < distribution.size(); i++)
		{
			distribution[i] /= sum;
		}
		return;
	}
	void to_csv(std::string fname, std::string seq = " "); //������ļ�
};

class Statistics
{
private:
	std::vector<double> data; //raw data
	double mean = 0; //��ֵ
	double median = 0; //��λ��
	double maxim = 0; //���ֵ
	double minimum = 0; //��Сֵ
	Distribution distribution; //�ֲ�

public:
	Statistics(std::vector<double> _data); //���캯��
	Distribution compute_distribution(double delta); //����ֲ�

	inline double get_mean() { return mean; } //����ƽ��ֵ
	inline double get_median() { return median; } //������λ��
	inline double get_maxim() { return maxim; } //�������ֵ
	inline double get_minimum() { return minimum; } //������Сֵ
	inline std::vector<double> get_data() { return data; } //����ԭʼ����
	inline Distribution get_distribution() { return distribution; } //���طֲ�
	inline void sort() //��С��������
	{
		std::sort(data.begin(), data.end());
	}
};
