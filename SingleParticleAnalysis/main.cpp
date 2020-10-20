#include "single_particle_analysis.h"
#include "configuration.h"
#include <iostream>
#include <string>
#include <fstream>


int main()
{
	std::cout << "enter equi data file name:";
	std::string fname;
	std::cin >> fname;
	Configuration config_preEqui(fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::pair);

	std::cout << "enter shear data file name: ";
	std::cin >> fname;
	Configuration_ParticleDynamic config_Equi(fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::single);
	config_Equi.compute_msd(config_preEqui);
	std::ofstream of_test;
	of_test.open("test.txt");
	for (size_t i = 0; i < config_Equi.get_msd().size(); i++)
	{
		of_test << config_Equi.get_msd()[i] << endl;
	}
	//std::cout << "enter pair sytle file name: ";
	//std::cin >> fname;
	//Configuration config_pair(fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::pair);
}