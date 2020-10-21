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
	minimum = *_data.begin();
	maxim = *_data.end();
	if (_data.size() % 2 == 0)
	{
		size_t size = _data.size();
		median = (_data[size / 2 - 1] + data[size / 2]) / 2;
	}
	else
	{
		median = _data[_data.size() / 2];
	}
	return;
}

Distribution Statistics::compute_distribution(size_t size)
{
	Distribution _dstbt;
	_dstbt.max = maxim;
	_dstbt.min = minimum;
	//_dstbt.size = size;
	double delta = (maxim - minimum) / size;
	_dstbt.distribution.resize(size, 0);
	for (size_t i = 0; i < data.size(); i++)
	{
		size_t index = floor((data[i] - minimum) / delta);
		_dstbt.distribution[index]++;
	}
}



void Distribution::to_file(std::string fname, std::string seq)
{
	std::ofstream ofile;
	ofile.open(fname);
	if (!ofile.is_open())
	{
		std::cerr << "File " << fname << " open failed!" << std::endl;
		throw ("File " + fname + " open failed");
	}
	ofile << "left boundary" << seq << "right boundarr" << seq << "distribution" << std::endl;
	for (size_t i = 0; i < distribution.size(); i++)
	{
		double delta = (max - min) / distribution.size();
		double left = i * delta;
		double right = i * delta + 1;
		ofile << left << seq << right << distribution[i];
		ofile << std::endl;
	}
	return;
}
