#include "statistics.h"

Statistics::Statistics(std::vector<double> _data)
{
	if (_data.empty())
	{
		std::cerr << "EMPTY data!" << std::endl;
		return;
	}
	data = _data;
	std::sort(_data.begin(), _data.end());
	minimum = _data.front();
	maxim = _data.back();

	if (_data.size() % 2 == 0)
	{
		double size = _data.size();
		median = (_data[size / 2 - 1] + data[size / 2]) / 2;
	}
	else
	{
		median = _data[_data.size() / 2];
	}
	mean = 0;
	for (size_t i = 0; i < _data.size(); i++)
	{
		mean += _data[i];
	}
	mean /= _data.size();
	return;
}

Distribution Statistics::compute_distribution(double delta)
{
	Distribution _dstbt;
	_dstbt.max = maxim;
	_dstbt.min = minimum;
	//double delta = (maxim - minimum) / (size - 1);
	size_t size = (maxim - minimum) / delta + 1;
	_dstbt.distribution.resize(size, 0);
	for (size_t i = 0; i < data.size(); i++)
	{
		size_t index = floor((data[i] - minimum) / delta); //将数据最小值对齐至0，再进行分布计算
		_dstbt.distribution[index]++;
	}
	distribution = _dstbt;
	return _dstbt;
}



void Distribution::to_csv(std::string fname, std::string seq)
{
	std::ofstream ofile;
	ofile.open(fname);
	if (!ofile.is_open())
	{
		//std::cerr << "File " << fname << " open failed!" << std::endl;
		throw ("File " + fname + " open failed");
	}
	ofile << "#" << "mid" << seq << "left boundary" << seq << "right boundarr" << seq << "distribution" << std::endl;
	for (size_t i = 0; i < distribution.size(); i++)
	{
		double delta = (max - min) / (distribution.size() - 1);
		double left = i * delta + min;
		double right = ((double)i + 1) * delta + min;
		double mid = (left + right) / 2;
		ofile << mid << seq << left << seq << right << seq << distribution[i];
		ofile << std::endl;
	}
	return;
}
