#include "single_particle_analysis.h"

void Configuration_StaticStructure::compute_CN(double r_cut)
{
	coordination_number.resize(get_particle().size());
	//盒子边界
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa_i = get_particle()[i];
		coordination_number[i].particle_id = p_pa_i.id;
		for (size_t j = 0; j < get_particle().size(); j++)
		{
			const Particle& p_pa_j = get_particle()[j];
			/*
			处理周期性坐标
			*/
			double dx = p_pa_i.rx - p_pa_j.rx;
			dx = (dx - std::floor(dx / lx * 2.) * lx / 2.);
			double dy = p_pa_i.ry - p_pa_j.ry;
			dy = (dy - std::floor(dy / ly * 2.) * ly / 2.);
			double dz = p_pa_i.rz - p_pa_j.rz;
			dz = (dz - std::floor(dz / lz * 2.) * lz / 2.);
			double dr = sqrt(dx * dx + dy * dy + dz * dz);
			if (dr < r_cut)
			{
				coordination_number[i].CN++;
			}
		}
	}
	for (size_t i = 0; i < coordination_number.size(); i++)
	{
		coordination_number[i].CN -= 1;
	}
	return;
}

const CN& Configuration_StaticStructure::get_CN(size_t _id) const
{
	for (size_t i = 0; i < coordination_number.size(); i++)
	{
		if (coordination_number[i].particle_id == _id)
		{
			return coordination_number[i];
		}
	}
	cerr << "particle " << _id << " CN not found!" << endl;
	throw ("particle " + std::to_string(_id) + " CN not found!");
}

void Configuration_StaticStructure::CN_to_file(std::string fname)
{
	vector<double> _vec_CN(get_particle_num());
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == coordination_number[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			_vec_CN[i] = coordination_number[i].CN;

		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			_vec_CN[i] = get_CN(get_particle()[i].id).CN;
		}

	}
	if (flag_same_seq)
	{
		to_dump(fname, { "CN" }, { _vec_CN });
	}
	else
	{
		string warning = "#WARNING: seq of vec_CN and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "CN" }, { _vec_CN }, { warning });
	}
	return;
}


void Configuration_ParticleDynamic::compute_msd(const Configuration& config_t0)
{
	//t时刻盒子边长
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();
	//0时刻盒子边长
	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_ylo();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zlo();
	//剪切情况下 l_t0 = l
	msd.resize(get_particle().size());
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa = get_particle()[i];
		const Particle& p_pa_t0 = config_t0.get_particle(p_pa.id);
		double drx = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
		double dry = p_pa.ry + p_pa.box_y * ly - p_pa_t0.ry - p_pa_t0.box_y * ly_t0;
		double drz = p_pa.rz + p_pa.box_z * lz - p_pa_t0.rz - p_pa_t0.box_z * lz_t0;
		double _msd = drx * drx + dry * dry + drz * drz;
		msd[i].particle_id = p_pa.id;
		msd[i].MSD_x = drx * drx;
		msd[i].MSD_y = dry * dry;
		msd[i].MSD_z = drz * drz;
		msd[i].MSD = _msd;
	}
	return;
}

const MSD& Configuration_ParticleDynamic::get_msd(size_t _id)
{

	for (size_t i = 0; i < msd.size(); i++)
	{
		if (msd[i].particle_id == _id)
		{
			return msd[i];
		}
	}
	cerr << "particle " << _id << " msd not found!" << endl;
	throw ("particle " + std::to_string(_id) + " msd not found!");
}

void Configuration_ParticleDynamic::compute_shear_MSDnonAffine(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, double step_time)
{

	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();

	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_ylo();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zlo();

	double dt = (get_timestep() - config_t0.get_timestep()) * step_time;

	msd_nonAffine.resize(get_particle().size() + 1);
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle p_pa = get_particle()[i];
		const Particle p_pa_t0 = config_t0.get_particle(p_pa.id);
		double drx = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
		double dry = p_pa.ry + p_pa.box_y * ly - p_pa_t0.ry - p_pa_t0.box_y * ly_t0;
		double drz = p_pa.rz + p_pa.box_z * lz - p_pa_t0.rz - p_pa_t0.box_z * lz_t0;

		switch (shear_direction)
		{
		case Configuration_ParticleDynamic::ShearDirection::xy:
			drx -= p_pa.ry * shear_rate * dt;
			break;
		case Configuration_ParticleDynamic::ShearDirection::xz:
			drx -= p_pa.rz * shear_rate * dt;
			break;
		case Configuration_ParticleDynamic::ShearDirection::yz:
			dry -= p_pa.rz * shear_rate * dt;
			break;
		default:
			break;
		}

		double _msd_nonAfine = drx * drx + dry * dry + drz * drz;
		msd_nonAffine[i].particle_id = p_pa.id;
		msd_nonAffine[i].MSD = _msd_nonAfine;
		msd_nonAffine[i].MSD_x = drx * drx;
		msd_nonAffine[i].MSD_y = dry * dry;
		msd_nonAffine[i].MSD_z = drz * drz;
	}
	return;
}

void Configuration_ParticleDynamic::compute_shear_MSDnonAffine(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, const vector<double>& gradient_ave_position, double step_time)
{
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();

	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_ylo();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zlo();

	double dt = (get_timestep() - config_t0.get_timestep()) * step_time;

	msd_nonAffine.resize(get_particle().size() + 1);
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle p_pa = get_particle()[i];
		const Particle p_pa_t0 = config_t0.get_particle(p_pa.id);
		double drx = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
		double dry = p_pa.ry + p_pa.box_y * ly - p_pa_t0.ry - p_pa_t0.box_y * ly_t0;
		double drz = p_pa.rz + p_pa.box_z * lz - p_pa_t0.rz - p_pa_t0.box_z * lz_t0;

		switch (shear_direction)
		{
		case Configuration_ParticleDynamic::ShearDirection::xy:
			drx -= gradient_ave_position[p_pa.id] * shear_rate * dt;
			break;
		case Configuration_ParticleDynamic::ShearDirection::xz:
			drx -= gradient_ave_position[p_pa.id] * shear_rate * dt;
			break;
		case Configuration_ParticleDynamic::ShearDirection::yz:
			dry -= gradient_ave_position[p_pa.id] * shear_rate * dt;
			break;
		default:
			break;
		}

		double _msd_nonAfine = drx * drx + dry * dry + drz * drz;
		msd_nonAffine[i].particle_id = p_pa.id;
		msd_nonAffine[i].MSD = _msd_nonAfine;
		msd_nonAffine[i].MSD_x = drx * drx;
		msd_nonAffine[i].MSD_y = dry * dry;
		msd_nonAffine[i].MSD_z = drz * drz;
	}
	return;
}

void Configuration_ParticleDynamic::compute_shear_MSDnonAffine_t0(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, double step_time)
{
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();

	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_ylo();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zlo();

	double dt = (get_timestep() - config_t0.get_timestep()) * step_time;

	msd_nonAffine.resize(get_particle().size() + 1);
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle p_pa = get_particle()[i];
		const Particle p_pa_t0 = config_t0.get_particle(p_pa.id);
		double drx = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
		double dry = p_pa.ry + p_pa.box_y * ly - p_pa_t0.ry - p_pa_t0.box_y * ly_t0;
		double drz = p_pa.rz + p_pa.box_z * lz - p_pa_t0.rz - p_pa_t0.box_z * lz_t0;

		switch (shear_direction)
		{
		case Configuration_ParticleDynamic::ShearDirection::xy:
			drx -= p_pa_t0.ry * shear_rate * dt;
			break;
		case Configuration_ParticleDynamic::ShearDirection::xz:
			drx -= p_pa_t0.rz * shear_rate * dt;
			break;
		case Configuration_ParticleDynamic::ShearDirection::yz:
			dry -= p_pa_t0.rz * shear_rate * dt;
			break;
		default:
			break;
		}

		double _msd_nonAfine = drx * drx + dry * dry + drz * drz;
		msd_nonAffine[i].particle_id = p_pa.id;
		msd_nonAffine[i].MSD = _msd_nonAfine;
		msd_nonAffine[i].MSD_x = drx * drx;
		msd_nonAffine[i].MSD_y = dry * dry;
		msd_nonAffine[i].MSD_z = drz * drz;
	}
	return;
}

const MSD& Configuration_ParticleDynamic::get_MSDnonAffine(size_t _id)
{
	for (size_t i = 0; i < msd.size(); i++)
	{
		if (msd_nonAffine[i].particle_id == _id)
		{
			return msd_nonAffine[i];
		}
	}
	cerr << "particle " << _id << " msd not found!" << endl;
	throw ("particle " + std::to_string(_id) + " msd not found!");
}

vector<Particle> Configuration_ParticleDynamic::pick_cross_gradient_boundary_particle(const Configuration& config_t0, ShearDirection shear_direction)
{
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa = get_particle()[i];
		const Particle& p_pa_t0 = config_t0.get_particle(p_pa.id);
		switch (shear_direction)
		{
		case Configuration_ParticleDynamic::ShearDirection::xy:
			if (p_pa.box_y != p_pa_t0.box_y)
			{
				cross_gradient_boundary_particle.push_back(p_pa);
			}
			break;
		case Configuration_ParticleDynamic::ShearDirection::xz:
			if (p_pa.box_z != p_pa_t0.box_z)
			{
				cross_gradient_boundary_particle.push_back(p_pa);
			}
			break;
		case Configuration_ParticleDynamic::ShearDirection::yz:
			if (p_pa.box_z != p_pa_t0.box_z)
			{
				cross_gradient_boundary_particle.push_back(p_pa);
			}
			break;
		default:
			break;
		}
	}
	return cross_gradient_boundary_particle;
}


void Configuration_ParticleDynamic::to_file_MSD(std::string fname)
{
	vector<double> _vec_msd(get_particle_num());
	//vector<double> _vec_msd_nonAffine(get_particle_num());
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == msd[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			_vec_msd[i] = msd[i].MSD;

		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			_vec_msd[i] = get_msd(get_particle()[i].id).MSD;
		}

	}
	if (flag_same_seq)
	{
		to_dump(fname, { "msd" }, { _vec_msd });
	}
	else
	{
		string warning = "#WARNING: seq of vec_msd and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "msd" }, { _vec_msd }, { warning });
	}
	return;
}

void Configuration_ParticleDynamic::to_file_nonAffineMSD(std::string fname)
{
	vector<double> _vec_msd_nonAffine(msd_nonAffine.size());
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == msd_nonAffine[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			_vec_msd_nonAffine[i] = msd_nonAffine[i].MSD;

		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			_vec_msd_nonAffine[i] = get_MSDnonAffine(get_particle()[i].id).MSD;
		}
	}
	if (flag_same_seq)
	{
		to_dump(fname, { "msd_nonAffine" }, { _vec_msd_nonAffine });
	}
	else
	{
		string warning = "#WARNING: seq of vec_msd_nonAffine and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "msd_nonAffine" }, { _vec_msd_nonAffine }, { warning });
	}
	return;
}

Configuration_ParticleDynamic Configuration_ParticleDynamic::gen_sub_config(const Configuration_ParticleDynamic& config_parents, vector<size_t> vec_id)
{
	Configuration_ParticleDynamic sub_config;
	sub_config = config_parents;
	sub_config.__clear_vec_particle();
	vector<Particle> vec_pa = config_parents.get_particle();
	for (size_t i = 0; i < vec_id.size(); i++)
	{
		sub_config.__add_particle(seek_id(vec_pa, vec_id[i]));
	}
	sub_config.set_particle_num(sub_config.get_particle().size());
	return sub_config;
}



