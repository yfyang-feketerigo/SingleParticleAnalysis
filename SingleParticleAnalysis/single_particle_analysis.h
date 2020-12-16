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
const double _STEP_TIME = 0.0025;
struct CN
{
	size_t particle_id;
	double cn;
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
	Configuration_StaticStructure()
		:Configuration() {};
	void compute_CN(double r_cut);
	inline const std::vector<CN>& get_CN()
	{
		return coordination_number;
	}
	const CN& get_CN(size_t _id) const;
	void CN_to_file(std::string fname);

};

struct MSD
{
	size_t particle_id;
	double MSD_x;
	double MSD_y;
	double MSD_z;
	double dx;
	double dy;
	double dz;
	double MSD;
};

struct Flow_deltaDisplacement
{
	size_t particle_id;
	double dr;
	int grad_box_change;
	size_t init_step;
	size_t now_step;
};

struct Flow_Ave_Velocity
{
	size_t particle_id;
	double ave_v;
	int grad_box_change;
};

class Configuration_ParticleDynamic :public Configuration_StaticStructure
{
	/*
	计算单粒子动态参量
	*/
private:
	std::vector<MSD> msd;
	std::vector<MSD> msd_nonAffine;
	std::vector<Particle> cross_gradient_boundary_particle;
	std::vector<Flow_deltaDisplacement> flow_delta_displacement;
	std::vector<Flow_Ave_Velocity> flow_ave_velocity;
public:
	/*Configuration_ParticleDynamic(std::string config_fname)
		:Configuration(config_fname) { };
		*/
		//Configuration_PrticleDynamic() {};
	Configuration_ParticleDynamic()
		:Configuration_StaticStructure() {};
	Configuration_ParticleDynamic(std::string config_fname, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single)
		:Configuration_StaticStructure(config_fname, _boxtype, _pairstyle) { };

	void compute_msd(const Configuration& config_t0);

	inline const std::vector<MSD>& get_MSD()
	{
		return msd;
	};
	const MSD& get_msd(size_t _id);

	enum class ShearDirection
	{
		xy,
		xz,
		yz
	};

	void compute_shear_MSDnonAffine(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, double step_time = _STEP_TIME); //按t时刻粒子梯度方向上的位置（周期性）计算仿射形变
	void compute_shear_MSDnonAffine(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, const vector<double>& gradient_ave_position, double step_time = _STEP_TIME);//按一段时间内梯度方向上的平均位置（周期性）计算仿射形变
	void compute_shear_MSDnonAffine_t0(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, double step_time = _STEP_TIME);
	void compute_shear_flow_displacement(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double step_time = _STEP_TIME);
	void compute_shear_flow_ave_velocity(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double step_time = _STEP_TIME);

	inline const std::vector<MSD>& get_MSDnonAffine()
	{
		return msd_nonAffine;
	};
	const MSD& get_MSDnonAffine(size_t _id);


	inline const std::vector<Flow_deltaDisplacement>& get_flow_displacement()
	{
		return flow_delta_displacement;
	}
	const Flow_deltaDisplacement& get_flow_displacement(size_t _id);

	inline const std::vector<Flow_Ave_Velocity>& get_flow_ave_velocity()
	{
		return flow_ave_velocity;
	}
	const Flow_Ave_Velocity& get_flow_ave_velocity(size_t _id);

	vector<Particle> pick_cross_gradient_boundary_particle(const Configuration& config_t0, ShearDirection shear_direction);

	void to_file_MSD(std::string fname);
	void to_file_nonAffineMSD(std::string fname);
	void to_file_flow_displacement(std::string fname);
	//void to_file_flow_displacement_non_cross_box(std::string fname);
	void to_file_flow_ave_velocity(std::string fname);
	//void to_file_flow_ave_velocity_non_cross_box(std::string fname);
	Configuration_ParticleDynamic gen_sub_config(const Configuration_ParticleDynamic& config_parents, vector<size_t> vec_id);
};