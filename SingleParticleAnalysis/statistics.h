/*

class Statistics :统计类，用于统计一维数据，可计算：
	最大值最小值，均值，中位数，分布

	struct Distribution：分布结构，储存分布数据类型。
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
	double min = 0; //数据中的最小值
	double max = 0; //数据中的最大值
	std::vector<double> distribution; //分布
	void normalization() //归一化
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
	void to_csv(std::string fname, std::string seq = " "); //输出到文件
};

class Statistics
{
private:
	std::vector<double> data; //raw data
	double mean = 0; //均值
	double median = 0; //中位数
	double maxim = 0; //最大值
	double minimum = 0; //最小值
	Distribution distribution; //分布

public:
	Statistics(std::vector<double> _data); //构造函数
	Distribution compute_distribution(double delta); //计算分布

	inline double get_mean() { return mean; } //返回平均值
	inline double get_median() { return median; } //返回中位数
	inline double get_maxim() { return maxim; } //返回最大值
	inline double get_minimum() { return minimum; } //返回最小值
	inline std::vector<double> get_data() { return data; } //返回原始数据
	inline Distribution get_distribution() { return distribution; } //返回分布
	inline void sort() //从小到大排序
	{
		std::sort(data.begin(), data.end());
	}
};
