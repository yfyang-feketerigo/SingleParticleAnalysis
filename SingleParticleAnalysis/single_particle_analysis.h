#pragma once
#include "input.h"
#include "particle.h"
#include "configuration.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
/*using std::string;
using std::cout;
using std::cerr;
using std::vector;
*/

struct CN
{
	size_t particle_id;
	double CN;
};
class Configuration_StaticStructure :public Configuration
{
	/*
	* 通过瞬态构象静态序参量
	* 继承Configuration类
	* 返回vector<double> 索引为粒子编号，值为配位数
	*/
private:
	//double z;
	std::vector<CN> coordination_number;
	//double
public:
	Configuration_StaticStructure(std::string config_fname)
		:Configuration(config_fname) { };

	Configuration_StaticStructure(std::string config_fname, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single)
		:Configuration(config_fname, _boxtype, _pairstyle) { };

	void compute_CN(double r_cut);
	inline std::vector<CN> get_CN()
	{
		return coordination_number;
	}
};

struct MSD
{
	size_t particle_id;
	double MSD_x;
	double MSD_y;
	double MSD_z;
	double MSD;
};

class Configuration_ParticleDynamic :public Configuration
{
	/*
	计算单粒子动态参量
	*/
private:
	std::vector<MSD> msd;
	std::vector<MSD> msd_nonAffine;
	std::vector<Particle> cross_gradient_boundary_particle;
public:
	Configuration_ParticleDynamic(std::string config_fname)
		:Configuration(config_fname) { };

	Configuration_ParticleDynamic(std::string config_fname, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single)
		:Configuration(config_fname, _boxtype, _pairstyle) { };

	void compute_msd(Configuration config_t0);

	inline std::vector<MSD> get_msd()
	{
		return msd;
	};

	enum class ShearDirection
	{
		xy,
		xz,
		yz
	};


	void compute_msd_nonAffine(Configuration config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, double dt = 0.0025);

	inline std::vector<MSD> get_msd_nonAffine()
	{
		return msd_nonAffine;
	};

	vector<Particle> pick_cross_gradient_boundary_particle(Configuration config_t0, ShearDirection shear_direction);
};