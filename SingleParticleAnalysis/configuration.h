//class, store configuration generated by lammps 'data' commend
#pragma once
#include "particle.h"
#include "input.h"
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <initializer_list>
//#include <algorithm>

using std::ofstream;
using std::clog;
using std::cout;
using std::to_string;
class Configuration
{
private:
	static const size_t LINE_SKIP_MAX = 2147483647;

	size_t HEAD_INFO_LINE; //data file description information, e.g. mass, pair etc.
	static size_t GAP_LINE;		 //description line between position and velocity

	size_t particle_num = 0; //total particle number
	unsigned long long timestep = 0;

	double xlo = 0;//xlo difined in lammps
	double ylo = 0;//ylo defined in lammps
	double zlo = 0;//zlo defined in lammps

	double xhi = 0;//xhi difined in lammps
	double yhi = 0;//yhi defined in lammps
	double zhi = 0;//zhi defined in lammps

	double xy = 0;//xy defined in lammps
	double xz = 0;//xz defined in lammps
	double yz = 0;//yz defined in lammps

	size_t type_num = 0;//number of particle types

	vector<string> strvec_mass_info; //store mass information in data file
	vector<string> strvec_pair_info; //store pair information in data file
	string str_atoms_info;
	vector<Particle> vec_particle; //particles container
	string filename;

public:
	enum class BoxType
	{
		orthogonal,
		tilt
	};
	enum class PairStyle
	{
		single,
		pair
	};
	Configuration(string config_file, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single);				   //

	inline size_t GET_LINE_MAX() { return LINE_SKIP_MAX; }
	inline size_t GET_HEAD_INFO_LINE() { return HEAD_INFO_LINE; }
	inline static size_t GET_GAP_LINE() { return GAP_LINE; }

	inline void SET_HEAD_INFO_LINE(size_t _HEAD_INFO_LINE) { HEAD_INFO_LINE = _HEAD_INFO_LINE; }
	inline static void SET_GAP_LINE(size_t _GAP_LINE) { GAP_LINE = _GAP_LINE; }

	inline void set_time_step(unsigned long long _timestep) { timestep = _timestep; };

	inline void set_particle_num(size_t _num) { particle_num = _num; } //set total number of particles
	inline size_t get_particle_num() const { return particle_num; }	   //return total number of particles
	inline size_t get_type_num() { return type_num; }				   //return number of particle types
	inline const vector<string>& get_mass_info() const				//return mass info
	{
		return strvec_mass_info;
	}
	inline const vector<string>& get_pair_info() const			//return pair info
	{
		return strvec_pair_info;
	}

	inline double get_xlo() const { return xlo; }//return xhi difined in lammps
	inline double get_xhi() const { return xhi; }//return xhi difined in lammps
	inline double get_ylo() const { return ylo; }//return ylo difined in lammps
	inline double get_yhi() const { return yhi; }//return yhi difined in lammps
	inline double get_zlo() const { return zlo; }//return zlo difined in lammps
	inline double get_zhi() const { return zhi; }//return zhi difined in lammps
	inline double get_xy() const { return xy; }//return xy difined in lammps
	inline double get_yz() const { return yz; }//return yz difined in lammps
	inline double get_xz() const { return xz; }//return xz difined in lammps

	inline double set_xboundary(double _xlo, double _xhi)
	{
		xlo = _xlo;
		xhi = _xhi;
	}
	inline double set_yboundary(double _ylo, double _yhi)
	{
		ylo = _ylo;
		yhi = _yhi;
	}
	inline double set_zboundary(double _zlo, double _zhi)
	{
		zlo = _zlo;
		zhi = _zhi;
	}
	inline double set_tilt(double _xy, double _xz, double _yz)
	{
		xy = _xy;
		xz = _xz;
		yz = _yz;
	}

	inline unsigned long long get_timestep() const { return timestep; }

	inline const string& get_filename() const { return filename; }

	inline const Particle& get_particle(size_t _id) const //return particle with given ID
	{
		for (size_t i = 0; i < vec_particle.size(); i++)
		{
			if (_id == vec_particle[i].id)
			{
				return vec_particle[i];
			}
		}
		std::cerr << "particle " << _id << " not found";
		throw("particle " + std::to_string(_id) + " not found");
	}
	inline const vector<Particle>& get_particle() const//return vector of all particles
	{
		return vec_particle;
	}


	void to_data(string fname, BoxType _boxtype = BoxType::tilt);
	void to_dump(string fname, std::initializer_list<string> add_para_name, std::initializer_list<vector<double>> add_para, vector<string> comments = {});//注意额外参量与粒子序号的对应关系
};

