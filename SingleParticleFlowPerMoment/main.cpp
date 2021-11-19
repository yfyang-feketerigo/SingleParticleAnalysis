
/*****************************************
20210406 add new pair style: none, for none pair info data file, modify settings to allow adjust pairsytle
of equi data file & shear data file
******************************************/
#include "../SingleParticleAnalysis/configuration.h"
#include "../SingleParticleAnalysis/single_particle_analysis.h"
#include <iostream>
#include <string>
#include <fstream>
#include <json/json.h>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

int main(int argc, char* argv[])
{
	try
	{
		bool flag_special_t0 = false;
		if (argc == 2 && string(argv[1]) == "help")
		{
			cout << "<./this_programm help> to see this message\n";
			cout << "excute like this: <./this_programm arg1 arg2 arg3[option]>\n";
			cout << "arg1: t0 time step\n";
			cout << "arg2: delta time step\n";
			cout << "arg3: optional, if set to \"flag_special_t0\", flowing parameters in jsonfile will be used:\n";
			cout << "    \"special_t0Equidata_boxtype\"\n";
			cout << "    \"special_t0Equidata_pairstyle\"\n";
		}
		else if (argc != 3 && argc != 4)
		{
			throw std::exception("illegal parameters!\n <. / this_programm help> to see help message\n");

		}
		else if (argc == 4)
		{
			string str_flag_special_t0(argv[3]);
			if ("flag_special_t0" == str_flag_special_t0) flag_special_t0 = true;
		}
		string str_timestep_t0(argv[1]);
		size_t timestep_t0 = std::stoll(str_timestep_t0);
		string str_delta_step(argv[2]);
		size_t delta_step = std::stoll(str_delta_step);




		std::cout << "timestep settings:\n";
		std::cout << "t0 timestep: " << timestep_t0 << "\n";
		std::cout << "delta timstep: " << delta_step << '\n';
		std::cout << "t timstep: " << delta_step + timestep_t0 << '\n';
		std::cout << "\n";
		std::ios_base::sync_with_stdio(false);
		std::cin.tie(NULL);
		boost::timer::auto_cpu_timer timer;
		Json::Value root;

		std::clog << "Reading PaAn_PerMoment.json..." << '\n';

		ifstream f_settings;
		f_settings.open("PaAn_PerMoment.json", std::ios_base::binary);
		if (!f_settings.is_open())
		{
			//std::cerr << "file \"PaAn_Settings.json\" open failed" << endl;
			throw std::exception("file \"PaAn_PerMoment.json\" open failed");
		}
		f_settings >> root;
		double wi = root["wi"].asDouble();
		double tau_alpha = root["tau_alpha"].asDouble();
		double rate = wi / tau_alpha;

		string fname_prefix = root["fname_prefix"].asString();
		string fname_postfix = root["fname_postfix"].asString();
		string data_fpath = root["data_fpath"].asString();
		string output_path = root["output_path"].asString();


		string str_t0data_pairstyle = root["t0data_pairstyle"].asString();
		if (flag_special_t0)
		{
			str_t0data_pairstyle = root["special_t0Equidata_pairstyle"].asString();
		}
		Configuration::PairStyle t0data_pairstyle = Configuration::PairStyle::none;
		if (str_t0data_pairstyle == "single")
		{
			t0data_pairstyle = Configuration::PairStyle::single;
		}
		else if (str_t0data_pairstyle == "pair")
		{
			t0data_pairstyle = Configuration::PairStyle::pair;
		}
		else if (str_t0data_pairstyle == "none")
		{
			t0data_pairstyle = Configuration::PairStyle::none;
		}
		else
		{
			throw std::exception(("WRONG pair style: " + str_t0data_pairstyle).c_str());
		}

		string str_t0data_boxtype = root["t0data_boxtype"].asString();
		if (flag_special_t0)
		{
			str_t0data_pairstyle = root["special_t0Equidata_boxtype"].asString();
		}
		auto t0data_boxsytle = Configuration::BoxType::orthogonal;
		if (str_t0data_boxtype == "orthogonal")
		{
			t0data_boxsytle = Configuration::BoxType::orthogonal;
		}
		else if (str_t0data_boxtype == "tilt")
		{
			t0data_boxsytle = Configuration::BoxType::tilt;
		}
		else
		{
			throw std::exception(("WRONG box type: " + str_t0data_boxtype).c_str());
		}


		string str_tdata_pairstyle = root["tdata_pairstyle"].asString();
		Configuration::PairStyle tdata_pairstyle = Configuration::PairStyle::none;
		if (str_tdata_pairstyle == "single")
		{
			tdata_pairstyle = Configuration::PairStyle::single;
		}
		else if (str_tdata_pairstyle == "pair")
		{
			tdata_pairstyle = Configuration::PairStyle::pair;
		}
		else if (str_tdata_pairstyle == "none")
		{
			tdata_pairstyle = Configuration::PairStyle::none;
		}
		else
		{
			throw std::exception(("WRONG pair style: " + str_tdata_pairstyle).c_str());
		}

		string str_tdata_boxtype = root["tdata_boxtype"].asString();
		auto tdata_boxsytle = Configuration::BoxType::orthogonal;
		if (str_tdata_boxtype == "orthogonal")
		{
			tdata_boxsytle = Configuration::BoxType::orthogonal;
		}
		else if (str_t0data_boxtype == "tilt")
		{
			tdata_boxsytle = Configuration::BoxType::tilt;
		}
		else
		{
			throw std::exception(("WRONG box type: " + str_t0data_boxtype).c_str());
		}




		boost::filesystem::path boost_path_check(data_fpath);
		if (!(boost::filesystem::exists(boost_path_check) && boost::filesystem::is_directory(boost_path_check)))
		{
			cerr << "data file path: " << boost_path_check << " not exits" << endl;
			throw std::exception(("data file path: " + boost_path_check.string() + " not exits").c_str());
		}

		void mkdir(std::string path);
		mkdir(output_path);

		bool flag_computeFlowDisplacement = root["computeFlowDisplacement"].asBool();
		bool flag_computeFlowAveVelocity = root["computeFlowAveVelocity"].asBool();


		string flow_displacement_dir = output_path + "FlowDisplacement/";
		string flow_displacement_non_cross_box_dir = output_path + "FlowDisplacement_nonCrossBox/";
		if (flag_computeFlowDisplacement)
		{
			mkdir(flow_displacement_dir);
		}

		string flow_ave_velocity_dir = output_path + "FlowAveVelocity/";
		string flow_ave_velocity_non_cross_box_dir = output_path + "FlowAveVelocity_nonCrossBox/";
		if (flag_computeFlowAveVelocity)
		{
			mkdir(flow_ave_velocity_dir);
		}

		string config_t0_fname = fname_prefix + to_string(timestep_t0) + fname_postfix;
		Configuration_StaticStructure config_t0
		(data_fpath + config_t0_fname, t0data_boxsytle, t0data_pairstyle);//Configuration::PairStyle::single);


		string str_timestep = to_string(timestep_t0 + delta_step);
		string config_t_fname = fname_prefix + str_timestep + fname_postfix;
		Configuration_ParticleDynamic config_t
		(data_fpath + config_t_fname, tdata_boxsytle, tdata_pairstyle);

		if (flag_computeFlowDisplacement)//!计算单粒子流场方向位移，以进一步计算流场问题，可以淘汰之前剪切带程序。
		{

			config_t.compute_shear_flow_displacement(config_t0, Configuration_ParticleDynamic::ShearDirection::xy);
			config_t.to_file_flow_displacement(flow_displacement_dir + "FlowDisplacement." + str_timestep);
		}

		if (flag_computeFlowAveVelocity)
		{

			config_t.compute_shear_flow_ave_velocity(config_t0, Configuration_ParticleDynamic::ShearDirection::xy);
			config_t.to_file_flow_ave_velocity(flow_ave_velocity_dir + "FlowAveVelocity." + str_timestep);
		}

		cout << "####################################################################" << '\n';
		return 0;
	}
	catch (const string e)
	{
		cout << e << '\n';
	}
	catch (const std::exception& e)
	{
		cout << "An exception occurred. Exception: " << e.what() << '\n';
	}
	catch (...)
	{
		cout << "An undifine exception occurred." << '\n';
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