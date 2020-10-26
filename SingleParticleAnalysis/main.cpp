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
		config_Equi.compute_CN(1.4);
		config_Equi.CN_to_file("test.CN.csv");
		ofstream of;

		std::unique_ptr<vector<double>> p_vec_temp(new vector<double>);
		p_vec_temp->resize(config_Equi.get_CN().size());
		for (size_t i = 0; i < p_vec_temp->size(); i++)
		{
			(*p_vec_temp)[i] = config_Equi.get_CN()[i].CN;
		}
		Statistics stat(*p_vec_temp);

		cout << stat.get_maxim() << " " << stat.get_minimum() << endl;
		cout << stat.get_mean() << endl;
		stat.compute_distribution(1);

		std::cout << "compute finished!" << std::endl;
		stat.get_distribution().to_csv("test_cn1.csv", ",");

		std::cout << "enter shear data file name: ";
		std::cin >> fname;
		Configuration_ParticleDynamic config_shear(fname, Configuration::BoxType::tilt);
		config_Equi.set_time_step(0);
		config_shear.compute_msd_nonAffine(config_Equi, Configuration_ParticleDynamic::ShearDirection::xy, (10 / 12000));
		of.close();

		config_shear.nonAffineMSD_to_file("test.nonAffineMSD.csv");
		Configuration_StaticStructure config_shear_transient(fname, Configuration::BoxType::tilt);
		config_shear_transient.compute_CN(0.94);
		config_shear_transient.CN_to_file("test.shear.CN.csv");

		return 0;
	}
	catch (const string e)
	{
		cout << e << endl;
	}

}

