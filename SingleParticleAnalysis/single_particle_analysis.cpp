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
		msd[i].dx = drx;
		msd[i].dy = dry;
		msd[i].dz = drz;
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

	msd_nonAffine.resize(get_particle().size());
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa = get_particle()[i];
		const Particle& p_pa_t0 = config_t0.get_particle(p_pa.id);
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
		msd_nonAffine[i].dx = drx;
		msd_nonAffine[i].dy = dry;
		msd_nonAffine[i].dz = drz;
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

	msd_nonAffine.resize(get_particle().size());
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa = get_particle()[i];
		const Particle& p_pa_t0 = config_t0.get_particle(p_pa.id);
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
		msd_nonAffine[i].dx = drx;
		msd_nonAffine[i].dy = dry;
		msd_nonAffine[i].dz = drz;
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

	msd_nonAffine.resize(get_particle().size());
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa = get_particle()[i];
		const Particle& p_pa_t0 = config_t0.get_particle(p_pa.id);
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
		msd_nonAffine[i].dx = drx;
		msd_nonAffine[i].dy = dry;
		msd_nonAffine[i].dz = drz;
	}
	return;
}

void Configuration_ParticleDynamic::compute_shear_flow_displacement(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double step_time)
{
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();

	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_ylo();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zlo();

	//double dt = (get_timestep() - config_t0.get_timestep()) * step_time;

	flow_delta_displacement.resize(get_particle_num() + 1);
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa = get_particle()[i];
		const Particle& p_pa_t0 = config_t0.get_particle(p_pa.id);
		flow_delta_displacement[i].particle_id = p_pa.id;
		flow_delta_displacement[i].init_step = config_t0.get_timestep();
		flow_delta_displacement[i].now_step = get_timestep();
		switch (shear_direction)
		{
		case Configuration_ParticleDynamic::ShearDirection::xy:
			flow_delta_displacement[i].dr = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
			flow_delta_displacement[i].grad_box_change = p_pa.box_y - p_pa_t0.box_y;
			break;
		case Configuration_ParticleDynamic::ShearDirection::xz:
			flow_delta_displacement[i].dr = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
			flow_delta_displacement[i].grad_box_change = p_pa.box_z - p_pa_t0.box_z;
			break;
		case Configuration_ParticleDynamic::ShearDirection::yz:
			flow_delta_displacement[i].dr = p_pa.ry + p_pa.box_y * ly - p_pa_t0.ry - p_pa_t0.box_y * ly_t0;
			flow_delta_displacement[i].grad_box_change = p_pa.box_z - p_pa_t0.box_z;
			break;
		default:
			break;
		}
	}
	return;
}

void Configuration_ParticleDynamic::compute_shear_flow_ave_velocity(const Configuration& config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double step_time)
{
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();

	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_ylo();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zlo();

	double dt = (get_timestep() - config_t0.get_timestep()) * step_time;

	flow_ave_velocity.resize(get_particle_num() + 1);
	for (size_t i = 0; i < get_particle().size(); i++)
	{
		const Particle& p_pa = get_particle()[i];
		const Particle& p_pa_t0 = config_t0.get_particle(p_pa.id);
		flow_ave_velocity[i].particle_id = p_pa.id;
		switch (shear_direction)
		{
		case Configuration_ParticleDynamic::ShearDirection::xy:
			flow_ave_velocity[i].ave_v = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
			flow_ave_velocity[i].grad_box_change = p_pa.box_y - p_pa_t0.box_y;
			break;
		case Configuration_ParticleDynamic::ShearDirection::xz:
			flow_ave_velocity[i].ave_v = p_pa.rx + p_pa.box_x * lx - p_pa_t0.rx - p_pa_t0.box_x * lx_t0;
			flow_ave_velocity[i].grad_box_change = p_pa.box_z - p_pa_t0.box_z;
			break;
		case Configuration_ParticleDynamic::ShearDirection::yz:
			flow_ave_velocity[i].ave_v = p_pa.ry + p_pa.box_y * ly - p_pa_t0.ry - p_pa_t0.box_y * ly_t0;
			flow_ave_velocity[i].grad_box_change = p_pa.box_z - p_pa_t0.box_z;
			break;
		default:
			break;
		}
		flow_ave_velocity[i].ave_v /= dt;
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

const Flow_deltaDisplacement& Configuration_ParticleDynamic::get_flow_displacement(size_t _id)
{
	for (size_t i = 0; i < flow_delta_displacement.size(); i++)
	{
		if (_id == flow_delta_displacement[i].particle_id)
		{
			return flow_delta_displacement[i];
		}
	}
	cerr << "particle " << _id << " flow displacement not found!" << endl;
	throw ("particle " + std::to_string(_id) + " flow displacement not found!");
}



const Flow_Ave_Velocity& Configuration_ParticleDynamic::get_flow_ave_velocity(size_t _id)
{
	for (size_t i = 0; i < msd.size(); i++)
	{
		if (_id == flow_ave_velocity[i].particle_id)
		{
			return flow_ave_velocity[i];
		}
	}
	cerr << "particle " << _id << " flow ave velocity not found!" << endl;
	throw ("particle " + std::to_string(_id) + " flow ave velocity not found!");
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
	vector<double> _vec_dx(get_particle_num());
	vector<double> _vec_dy(get_particle_num());
	vector<double> _vec_dz(get_particle_num());
	//vector<double> _vec_msd_nonAffine(get_particle_num());
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == msd[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			_vec_msd[i] = msd[i].MSD;
			_vec_dx[i] = msd[i].dx;
			_vec_dy[i] = msd[i].dy;
			_vec_dz[i] = msd[i].dz;

		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			_vec_msd[i] = get_msd(get_particle()[i].id).MSD;
			_vec_dx[i] = msd[i].dx;
			_vec_dy[i] = msd[i].dy;
			_vec_dz[i] = msd[i].dz;
		}

	}
	if (flag_same_seq)
	{
		to_dump(fname, { "msd","dx","dy","dz" }, { _vec_msd,_vec_dx,_vec_dy,_vec_dz });
	}
	else
	{
		string warning = "#WARNING: seq of vec_msd and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "msd","dx","dy","dz" }, { _vec_msd,_vec_dx,_vec_dy,_vec_dz });
	}
	return;
}

void Configuration_ParticleDynamic::to_file_nonAffineMSD(std::string fname)
{
	vector<double> _vec_msd_nonAffine(get_particle().size());
	vector<double> _vec_nonAffine_dx(get_particle().size());
	vector<double> _vec_nonAffine_dy(get_particle().size());
	vector<double> _vec_nonAffine_dz(get_particle().size());
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == msd_nonAffine[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			_vec_msd_nonAffine[i] = msd_nonAffine[i].MSD;
			_vec_nonAffine_dx[i] = msd_nonAffine[i].dx;
			_vec_nonAffine_dy[i] = msd_nonAffine[i].dy;
			_vec_nonAffine_dz[i] = msd_nonAffine[i].dz;

		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			_vec_msd_nonAffine[i] = get_MSDnonAffine(get_particle()[i].id).MSD;
			_vec_nonAffine_dx[i] = msd_nonAffine[i].dx;
			_vec_nonAffine_dy[i] = msd_nonAffine[i].dy;
			_vec_nonAffine_dz[i] = msd_nonAffine[i].dz;
		}
	}
	if (flag_same_seq)
	{
		to_dump(fname, { "msd_nonAffine" }, { _vec_msd_nonAffine, _vec_nonAffine_dx,_vec_nonAffine_dy,_vec_nonAffine_dz });
	}
	else
	{
		string warning = "#WARNING: seq of vec_msd_nonAffine and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "msd_nonAffine" "dx","dy","dz" }, { _vec_msd_nonAffine, _vec_nonAffine_dx,_vec_nonAffine_dy,_vec_nonAffine_dz }, { warning });
	}
	return;
}

void Configuration_ParticleDynamic::to_file_flow_displacement(std::string fname)
{
	//todo 尚未测试
	vector<double> _vec_flow_displacement(flow_delta_displacement.size());
	vector<double> _vec_flow_grad_box_change(flow_delta_displacement.size());
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == flow_delta_displacement[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			_vec_flow_displacement[i] = flow_delta_displacement[i].dr;
			_vec_flow_grad_box_change[i] = flow_delta_displacement[i].grad_box_change;
		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			_vec_flow_displacement[i] = get_flow_displacement(get_particle()[i].id).dr;
			_vec_flow_grad_box_change[i] = get_flow_displacement(get_particle()[i].id).grad_box_change;
		}
	}
	string step_info = "# init_step: " + to_string(flow_delta_displacement[0].init_step);
	step_info += " now_step: " + to_string(flow_delta_displacement[0].now_step);
	if (flag_same_seq)
	{
		to_dump(fname, { "flow_delta_displacement", "grad_box_change" }, { _vec_flow_displacement, _vec_flow_grad_box_change }, { step_info });
	}
	else
	{
		string warning = "#WARNING: seq of vec_flow_displacement and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "flow_delta_displacement", "grad_box_change" }, { _vec_flow_displacement, _vec_flow_grad_box_change }, { warning,step_info });
	}
}
/*
void Configuration_ParticleDynamic::to_file_flow_displacement_non_cross_box(std::string fname)
{
	vector<double> _vec_flow_displacement;
	vector<double> _vec_flow_grad_box_change;
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == flow_displacement[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			if (flow_displacement[i].grad_box_change == 0)
			{
				_vec_flow_displacement.push_back(flow_displacement[i].dr);
				_vec_flow_grad_box_change.push_back(flow_displacement[i].grad_box_change);
			}
		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			if (get_flow_displacement(get_particle()[i].id).grad_box_change == 0)
			{
				_vec_flow_displacement.push_back(get_flow_displacement(get_particle()[i].id).dr);
				_vec_flow_grad_box_change.push_back(get_flow_displacement(get_particle()[i].id).grad_box_change);
			}
		}
	}
	if (flag_same_seq)
	{
		to_dump(fname, { "flow_displacement", "grad_box_change" }, { _vec_flow_displacement, _vec_flow_grad_box_change });
	}
	else
	{
		string warning = "#WARNING: seq of vec_flow_displacement and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "flow_displacement", "grad_box_change" }, { _vec_flow_displacement, _vec_flow_grad_box_change }, { warning });
	}
}
*/
void Configuration_ParticleDynamic::to_file_flow_ave_velocity(std::string fname)
{
	//todo 尚未测试
	vector<double> _vec_flow_ave_velocity(flow_ave_velocity.size());
	vector<double> _vec_flow_grad_box_change(flow_ave_velocity.size());
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == flow_ave_velocity[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			_vec_flow_ave_velocity[i] = flow_ave_velocity[i].ave_v;
			_vec_flow_grad_box_change[i] = flow_ave_velocity[i].grad_box_change;
		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			_vec_flow_ave_velocity[i] = get_flow_ave_velocity(get_particle()[i].id).ave_v;
			_vec_flow_grad_box_change[i] = get_flow_ave_velocity(get_particle()[i].id).grad_box_change;
		}
	}
	if (flag_same_seq)
	{
		to_dump(fname, { "flow_delta_displacement", "grad_box_change" }, { _vec_flow_ave_velocity, _vec_flow_grad_box_change });
	}
	else
	{
		string warning = "#WARNING: seq of vec_flow_ave_veclocity and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "flow_ave_velocity" , "grad_box_change" }, { _vec_flow_ave_velocity, _vec_flow_grad_box_change }, { warning });
	}
}

/*void Configuration_ParticleDynamic::to_file_flow_ave_velocity_non_cross_box(std::string fname)
{
	vector<double> _vec_flow_ave_velocity;
	vector<double> _vec_flow_grad_box_change;
	bool flag_same_seq = true;

	for (size_t i = 0; i < get_particle().size(); i++)
	{
		if (get_particle()[i].id == flow_ave_velocity[i].particle_id)
		{
			flag_same_seq = flag_same_seq && true;
			if (flow_ave_velocity[i].grad_box_change == 0)
			{
				_vec_flow_ave_velocity.push_back(flow_ave_velocity[i].ave_v);
				_vec_flow_grad_box_change.push_back(flow_ave_velocity[i].grad_box_change);
			}
		}
		else
		{
			flag_same_seq = flag_same_seq && false;
			if (get_flow_ave_velocity(get_particle()[i].id).grad_box_change == 0)
			{
				_vec_flow_ave_velocity.push_back(get_flow_ave_velocity(get_particle()[i].id).ave_v);
				_vec_flow_grad_box_change.push_back(get_flow_ave_velocity(get_particle()[i].id).grad_box_change);
			}
		}
	}
	if (flag_same_seq)
	{
		to_dump(fname, { "flow_displacement", "grad_box_change" }, { _vec_flow_ave_velocity, _vec_flow_grad_box_change });
	}
	else
	{
		string warning = "#WARNING: seq of vec_flow_ave_veclocity and vec_particle not match";
		std::clog << warning << endl;
		to_dump(fname, { "flow_ave_velocity" , "grad_box_change" }, { _vec_flow_ave_velocity, _vec_flow_grad_box_change }, { warning });
	}
}*/


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



