#include "single_particle_analysis.h"

void Configuration_StaticStructure::compute_CN(double r_cut)
{
	coordination_number.resize(get_particle().size());
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();
	for (size_t i = 0; i < get_particle().size() - 1; i++)
	{
		//const std::unique_ptr<Particle> p_pa_i(&get_particle()[i]);
		const Particle& p_pa_i = get_particle()[i];
		coordination_number[i].particle_id = p_pa_i.id;
		for (size_t j = i + 1; j < get_particle().size(); j++)
		{
			const Particle& p_pa_j = get_particle()[j];
			coordination_number[j].particle_id = p_pa_j.id;
			double dx = p_pa_i.rx - p_pa_j.rx;
			dx = dx - std::floor(dx / lx * 2) * lx / 2;
			double dy = p_pa_i.ry - p_pa_j.ry;
			dy = dy - std::floor(dy / ly * 2) * ly / 2;
			double dz = p_pa_i.rz - p_pa_j.rz;
			dz = dz - std::floor(dz / lz * 2) * lz / 2;
			double dr = sqrt(dx * dx + dy * dy + dz * dz);
			if (dr < r_cut)
			{
				coordination_number[i].CN++;
				coordination_number[j].CN++;
			}
		}
	}
	return;
}

void Configuration_ParticleDynamic::compute_msd(Configuration config_t0)
{

	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();

	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_yhi();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zhi();

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

void Configuration_ParticleDynamic::compute_msd_nonAffine(Configuration config_t0, Configuration_ParticleDynamic::ShearDirection shear_direction, double shear_rate, double step_time)
{

	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();

	double lx_t0 = config_t0.get_xhi() - config_t0.get_xlo();
	double ly_t0 = config_t0.get_yhi() - config_t0.get_yhi();
	double lz_t0 = config_t0.get_zhi() - config_t0.get_zhi();

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

vector<Particle> Configuration_ParticleDynamic::pick_cross_gradient_boundary_particle(Configuration config_t0, ShearDirection shear_direction)
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


