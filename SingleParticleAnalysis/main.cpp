#include "single_particle_analysis.h"
#include "configuration.h"
#include <iostream>
#include <string>
#include <fstream>


int main()
{
	//std::cout << "enter preEqui data file name:";
	std::string fname;
	//std::cin >> fname;
	//Configuration config_tilt(fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::pair);

	std::cout << "enter equi data file name: ";
	std::cin >> fname;
	Configuration_ParticleDynamic config_Equi(fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::single);
	config_Equi.set_time_step(0);
	std::cout << "enter shear data file name: ";
	std::cin >> fname;
	Configuration_ParticleDynamic config_shear(fname, Configuration::BoxType::tilt, Configuration::PairStyle::single);
	double shear_rate = 0.1 / 12000;
	config_shear.compute_msd_nonAffine(config_Equi, Configuration_ParticleDynamic::ShearDirection::xy, shear_rate);

	std::ofstream of_test;
	of_test.open("test.txt");

	double counter = 0;
	for (size_t i = 0; i < config_shear.get_particle().size(); i++)
	{
		const Particle* p_pa = &config_shear.get_particle()[i];
		const Particle* p_pa_t0 = &config_Equi.get_particle(p_pa->id);
		if (p_pa->box_y != p_pa_t0->box_y)
		{
			counter++;
			cout << p_pa->id << " " << p_pa->box_y << " " << p_pa_t0->box_y << endl;;
		}
	}
	cout << "###################" << endl;
	cout << counter << endl;
	cout << "###################" << endl;
	for (size_t i = 0; i < config_shear.get_msd_nonAffine().size(); i++)
	{
		of_test << config_shear.get_msd_nonAffine()[i] << endl;
	}
	//std::cout << "enter pair sytle file name: ";
	//std::cin >> fname;
	//Configuration config_pair(fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::pair);
}