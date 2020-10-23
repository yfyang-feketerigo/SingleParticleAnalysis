#include "single_particle_analysis.h"
#include "configuration.h"
#include "statistics.h"
#include <iostream>
#include <string>
#include <fstream>


int main()
{
	try
	{
		std::string fname;
		std::cout << "enter equi data file name: ";
		std::cin >> fname;
		Configuration_StaticStructure config_Equi(fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::single);
		config_Equi.set_time_step(0);
		config_Equi.compute_CN(0.94);
		ofstream of;
		of.open("test_CN.txt");
		for (size_t i = 0; i < config_Equi.get_CN().size(); i++)
		{
			of << config_Equi.get_CN()[i].particle_id << " " << config_Equi.get_CN()[i].CN << endl;
		}


		std::unique_ptr<vector<double>> p_vec_temp(new vector<double>);
		p_vec_temp->resize(config_Equi.get_CN().size());
		for (size_t i = 0; i < p_vec_temp->size(); i++)
		{
			(*p_vec_temp)[i] = config_Equi.get_CN()[i].CN;
		}
		Statistics stat(*p_vec_temp);
		//p_vec_temp.reset();
		//p_vec_temp.get_deleter();

		cout << stat.get_maxim() << " " << stat.get_minimum() << endl;
		cout << stat.get_mean() << endl;


		stat.compute_distribution(1);

		std::cout << "compute finished!" << std::endl;
		stat.get_distribution().to_csv("test_cn1.csv", ",");

		std::cout << "enter shear data file name: ";
		std::cin >> fname;
		Configuration_ParticleDynamic config_shear(fname, Configuration::BoxType::tilt);
		config_shear.compute_msd_nonAffine(config_Equi, Configuration_ParticleDynamic::ShearDirection::xy, (10 / 12000));
		config_Equi.set_time_step(0);
		of.close();
		of.open("test_msd.csv");
		for (size_t i = 0; i < config_shear.get_msd_nonAffine().size(); i++)
		{
			if (config_shear.get_msd_nonAffine()[i].MSD > 400)
			{
				size_t _id = config_shear.get_msd_nonAffine()[i].particle_id;
				of << _id << " ";
				const Particle& pa = config_shear.get_particle(_id);
				of << pa.type << " " << pa.rx << " " << pa.ry << " " << pa.rz;
				of << endl;
			}
		}

		return 0;
	}
	catch (const string e)
	{
		cout << e << endl;
	}

}

