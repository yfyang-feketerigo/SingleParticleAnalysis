#include "single_particle_analysis.h"
#include "configuration.h"
//#include "statistics.h"
#include <iostream>
#include <string>
#include <fstream>
#include <json/json.h>
//#include <boost/timer.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

int main()
{
	try
	{
		boost::timer::auto_cpu_timer timer;
		Json::Value root;

		std::clog << "Reading PaAn_Settings.json..." << endl;

		ifstream f_settings;
		f_settings.open("PaAn_Settings.json", std::ios_base::binary);
		if (!f_settings.is_open())
		{
			std::cerr << "file \"PaAn_Settings.json\" open failed" << endl;
			throw "file \"PaAn_Settings.json\" open failed";
		}
		f_settings >> root;
		double wi = root["wi"].asDouble();
		double tau_alpha = root["tau_alpha"].asDouble();
		double rate = wi / tau_alpha;
		string equi_fname = root["equi_config_fname"].asString();
		//size_t start_step = root["start_step"].asUInt();
		//size_t end_step = root["end_step"].asLargestUInt();
		size_t start_moment = root["start_moment"].asLargestUInt();
		size_t moment_number = root["moment_number"].asLargestUInt();
		size_t delta_step = root["delta_step"].asLargestUInt();
		string fname_prefix = root["fname_prefix"].asString();
		string fname_postfix = root["fname_postfix"].asString();
		string data_fpath = root["data_fpath"].asString();
		string output_path = root["output_path"].asString();
		string equi_data_fpath = root["equi_data_fpath"].asString();
		double CN_rcut = root["CN_rcut"].asDouble();

		boost::filesystem::path boost_path_check(data_fpath);
		if (!(boost::filesystem::exists(boost_path_check) && boost::filesystem::is_directory(boost_path_check)))
		{
			cerr << "data file path: " << boost_path_check << " not exits" << endl;
			throw ("data file path: " + boost_path_check.string() + " not exits");
		}
		boost_path_check = boost::filesystem::path(equi_data_fpath);
		if (!(boost::filesystem::exists(boost_path_check) && boost::filesystem::is_directory(boost_path_check)))
		{
			cerr << "data file path: " << boost_path_check << " not exits" << endl;
			throw ("data file path: " + boost_path_check.string() + " not exits");
		}

		void mkdir(std::string path);
		mkdir(output_path);

		//string statisticDir = output_path + "stastics/";

		Configuration_StaticStructure config_equi
		(equi_data_fpath + equi_fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::single);
		config_equi.set_time_step(0); //平衡态设为0时间点

		bool flag_CN = root["computeCN"].asBool();
		bool flag_MSD = root["computeMSD"].asBool();
		bool flag_MSDnonAffine = root["computeMSDnonAffine"].asBool();
		bool flag_MSDnonAffine_t0 = root["MSDnonAffine_t0"].asBool();
		bool flag_MSDnonAffine_ave_gradient = root["MSDnonAffine_ave_gradient"].asBool();

		string CNdir = output_path + "CN/";
		if (flag_CN) mkdir(CNdir);

		string MSDdir = output_path + "MSD/";
		if (flag_MSD) mkdir(MSDdir);

		string MSDnonAffinedir = output_path + "MSDnonAffine/";
		if (flag_MSDnonAffine) mkdir(MSDnonAffinedir);

		string MSDnonAffine_gradient_ave_dir = output_path + "MSDnonAffineGrdAve/";
		if (flag_MSDnonAffine_ave_gradient) mkdir(MSDnonAffine_gradient_ave_dir);

		string MSDnonAffine_t0_dir = output_path + "MSDnonAffine_t0/";
		if (flag_MSDnonAffine_t0) mkdir(MSDnonAffine_t0_dir);

		if (flag_CN)
		{
			config_equi.compute_CN(CN_rcut);
			config_equi.CN_to_file(CNdir + "CN.0");
		}
		vector<double> y_sum(config_equi.get_particle().size(), 0);
		vector<double> x_sum(config_equi.get_particle().size(), 0);
		vector<double> z_sum(config_equi.get_particle().size(), 0);
		for (size_t i = 0; i < config_equi.get_particle().size(); i++)
		{
			const Particle& pa = config_equi.get_particle()[i];
			x_sum[pa.id] = pa.rx;
			y_sum[pa.id] = pa.ry;
			z_sum[pa.id] = pa.rz;
		}

		size_t ncounter = 0;
		for (size_t imoment = start_moment; imoment <= moment_number; imoment++)
		{
			size_t istep = imoment * delta_step;
			ncounter++;
			string str_istep = to_string(istep);
			string config_t_fname = data_fpath + fname_prefix + str_istep + fname_postfix;
			Configuration_ParticleDynamic config_t(config_t_fname, Configuration::BoxType::tilt);
			if (flag_CN)
			{
				config_t.compute_CN(CN_rcut);
				config_t.CN_to_file(CNdir + "CN." + str_istep);
			}

			if (flag_MSD)
			{
				config_t.compute_msd(config_equi);
				config_t.to_file_MSD(MSDdir + "MSD." + str_istep);
			}

			if (flag_MSDnonAffine)
			{
				config_t.compute_shear_MSDnonAffine(config_equi, Configuration_ParticleDynamic::ShearDirection::xy, rate);
				config_t.to_file_nonAffineMSD(MSDnonAffinedir + "MSDnonAffine." + str_istep);
			}

			if (flag_MSDnonAffine_t0)
			{
				config_t.compute_shear_MSDnonAffine_t0(config_equi, Configuration_ParticleDynamic::ShearDirection::xy, rate);
				config_t.to_file_nonAffineMSD(MSDnonAffine_t0_dir + "MSDnonAffine_t0." + str_istep);
			}

			if (flag_MSDnonAffine_ave_gradient)
			{
				for (size_t i = 0; i < config_t.get_particle().size(); i++)
				{
					const Particle& pa = config_t.get_particle()[i];
					x_sum[pa.id] += pa.rx;
					y_sum[pa.id] += pa.ry;
					z_sum[pa.id] += pa.rz;
				}
				vector<double> y_ave = y_sum;
				for (size_t i = 0; i < y_ave.size(); i++)
				{
					y_ave[i] /= ncounter;
				}
				config_t.compute_shear_MSDnonAffine(config_equi, Configuration_ParticleDynamic::ShearDirection::xy, rate, y_ave);
				config_t.to_file_nonAffineMSD(MSDnonAffine_gradient_ave_dir + "MSDnonAffineGrdAve." + str_istep);
			}

		}
		cout << "####################################################################" << endl;
		//cout << "PROGRAM END, TOTAL TIME USED: " << timer.elapsed() << "s" << endl;
		return 0;
	}
	catch (const string e)
	{
		cout << e << endl;
	}

}

void mkdir(std::string path)
{
	boost::filesystem::path bpath(path);
	if (!boost::filesystem::exists(path))
	{
		boost::filesystem::create_directories(bpath);
	}
}