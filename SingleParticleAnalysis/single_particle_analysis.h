#pragma once
#include "input.h"
#include "particle.h"
#include "configuration.h"
#include <string>
#include <iostream>
#include <vector>
/*using std::string;
using std::cout;
using std::cerr;
using std::vector;
*/

class Configuration_StaticStructure :public Configuration
{
	/*
	* 通过瞬态构象静态序参量
	* 继承Configuration类
	* 返回vector<double> 索引为粒子编号，值为配位数
	*/
private:
	//double z;
	std::vector<unsigned int> coordination_number;
	//double
public:
	Configuration_StaticStructure(std::string config_fname)
		:Configuration(config_fname) { };

	Configuration_StaticStructure(std::string config_fname, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single)
		:Configuration(config_fname, _boxtype, _pairstyle) { };

	void compute_CN(double r_cut);
	inline std::vector<unsigned int> get_CN()
	{
		return coordination_number;
	}
};

class Configuration_ParticleDynamic :public Configuration
{
private:
	std::vector<double> msd;
	std::vector<double> msd_nonAffine;
public:
	Configuration_ParticleDynamic(std::string config_fname)
		:Configuration(config_fname) { };

	Configuration_ParticleDynamic(std::string config_fname, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single)
		:Configuration(config_fname, _boxtype, _pairstyle) { };

	void compute_msd(Configuration config_t0);

	inline std::vector<double> get_msd()
	{
		return msd;
	};

	enum class ShearDirection
	{
		xy,
		xz,
		yz
	};

	void compute_msd_nonAffine(Configuration config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction);

	inline std::vector<double> get_msd_nonAffine()
	{
		return msd_nonAffine;
	};

};