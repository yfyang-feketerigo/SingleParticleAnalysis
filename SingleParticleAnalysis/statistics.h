#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
struct Distribution
{
	double min = 0;
	double max = 0;
	std::vector<double> distribution;
	void normalization()
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
	void to_file(std::string fname, std::string seq = " ");
};

class Statistics
{
private:
	std::vector<double> data;
	double mean = 0;
	double median = 0;
	double maxim = 0;
	double minimum = 0;
	Distribution distribution;

public:
	Statistics(std::vector<double> _data);
	Distribution compute_distribution(size_t size);

	inline double get_mean() { return mean; }
	inline double get_median() { return median; }
	inline double get_max() { return maxim; }
	inline double get_minimum() { return minimum; }
	inline std::vector<double> get_data() { return data; }
	inline Distribution get_distribution() { return distribution; }
	inline void sort() {
		std::sort(data.begin(), data.end());
	}
};

