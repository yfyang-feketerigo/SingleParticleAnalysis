
/*****************************************
20210406 add new pair style: none, for none pair info data file, modify settings to allow adjust pairsytle
of equi data file & shear data file
******************************************/
#include "single_particle_analysis.h"
#include "configuration.h"
#include <iostream>
#include <string>
#include <fstream>
#include <json/json.h>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

int main()
{
	using namespace std;
	using namespace Input;
	try
	{
		std::ios_base::sync_with_stdio(false);
		std::cin.tie(NULL);
		boost::timer::auto_cpu_timer timer;
		Json::Value root;

		std::clog << "Reading PaAn_Settings.json..." << '\n';

		ifstream f_settings;
		f_settings.open("PaAn_Settings.json", std::ios_base::binary);
		if (!f_settings.is_open())
		{
			//std::cerr << "file \"PaAn_Settings.json\" open failed" << endl;
			throw std::runtime_error("file \"PaAn_Settings.json\" open failed");
		}
		f_settings >> root;
		string str_shear_direction = root["shear_direction"].asString();
		auto shear_direction = Configuration_ParticleDynamic::ShearDirection::xy;
		if ("xy" == str_shear_direction)
		{
			shear_direction = Configuration_ParticleDynamic::ShearDirection::xy;
		}
		else if ("yz" == str_shear_direction)
		{
			shear_direction = Configuration_ParticleDynamic::ShearDirection::yz;
		}
		else if ("xz" == str_shear_direction)
		{
			shear_direction = Configuration_ParticleDynamic::ShearDirection::xz;
		}
		else
		{
			string err_msg = "illegal shear direction: " + str_shear_direction;
			throw std::runtime_error(err_msg.c_str());
		}

		double wi = root["wi"].asDouble();
		double tau_alpha = root["tau_alpha"].asDouble();
		double rate = wi / tau_alpha;
		string t0_fname = root["t0config_fname"].asString();
		size_t start_moment = root["start_moment"].asLargestUInt();
		size_t moment_number = root["moment_number"].asLargestUInt();
		size_t delta_step = root["delta_step"].asLargestUInt();
		size_t t0_shift_timestep = root["t0_shift_step"].asLargestUInt();
		string fname_prefix = root["fname_prefix"].asString();
		string fname_postfix = root["fname_postfix"].asString();
		string data_fpath = root["data_fpath"].asString();
		string output_path = root["output_path"].asString();
		string t0data_fpath = root["t0data_fpath"].asString();

		string str_t0data_pairstyle = root["t0data_pairstyle"].asString();
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
			//cerr << ("WRONG pair style: " + str_t0data_pairstyle) << endl;
			throw std::runtime_error("WRONG pair style: " + str_t0data_pairstyle);
		}

		string str_t0data_boxtype = root["t0data_boxtype"].asString();
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
			cerr << ("WRONG box type: " + str_t0data_boxtype) << endl;
			throw std::runtime_error("WRONG box type: " + str_t0data_boxtype);
		}

		string str_shear_data_pairstyle = root["shear_data_pairstyle"].asString();
		Configuration::PairStyle shear_data_pairstyle = Configuration::PairStyle::none;
		if (str_shear_data_pairstyle == "single")
		{
			shear_data_pairstyle = Configuration::PairStyle::single;
		}
		else if (str_shear_data_pairstyle == "pair")
		{
			shear_data_pairstyle = Configuration::PairStyle::pair;
		}
		else if (str_shear_data_pairstyle == "none")
		{
			shear_data_pairstyle = Configuration::PairStyle::none;
		}
		else
		{
			cerr << ("WRONG pair style: " + str_shear_data_pairstyle) << endl;
			throw std::runtime_error("WRONG pair style: " + str_shear_data_pairstyle);
		}


		double CN_rcut = root["CN_rcut"].asDouble();


		boost::filesystem::path boost_path_check(data_fpath);
		if (!(boost::filesystem::exists(boost_path_check) && boost::filesystem::is_directory(boost_path_check)))
		{
			cerr << "data file path: " << boost_path_check << " not exits" << endl;
			throw std::runtime_error("data file path: " + boost_path_check.string() + " not exits");
		}
		boost_path_check = boost::filesystem::path(t0data_fpath);
		if (!(boost::filesystem::exists(boost_path_check) && boost::filesystem::is_directory(boost_path_check)))
		{
			cerr << "data file path: " << boost_path_check << " not exits" << endl;
			throw std::runtime_error("data file path: " + boost_path_check.string() + " not exits");
		}

		void mkdir(std::string path);
		mkdir(output_path);


		Configuration_StaticStructure config_t0
		(t0data_fpath + t0_fname, t0data_boxsytle, t0data_pairstyle);//Configuration::PairStyle::single);
		config_t0.set_time_step(0 + t0_shift_timestep); //设定t0

		bool flag_CN = root["computeCN"].asBool();
		bool flag_MSD = root["computeMSD"].asBool();
		bool flag_MSDnonAffine = root["computeMSDnonAffine"].asBool();
		bool flag_MSDnonAffine_t0 = root["MSDnonAffine_t0"].asBool();
		bool flag_MSDnonAffine_ave_gradient = root["MSDnonAffine_ave_gradient"].asBool();
		bool flag_computeFlowDisplacement = root["computeFlowDisplacement"].asBool();
		bool flag_computeFlowAveVelocity = root["computeFlowAveVelocity"].asBool();
		string CNdir = output_path + "cn/";
		if (flag_CN) mkdir(CNdir);

		string MSDdir = output_path + "MSD/";
		if (flag_MSD) mkdir(MSDdir);

		string MSDnonAffinedir = output_path + "MSDnonAffine/";
		if (flag_MSDnonAffine) mkdir(MSDnonAffinedir);

		string MSDnonAffine_gradient_ave_dir = output_path + "MSDnonAffineGrdAve/";
		if (flag_MSDnonAffine_ave_gradient) mkdir(MSDnonAffine_gradient_ave_dir);

		string MSDnonAffine_t0_dir = output_path + "MSDnonAffine_t0/";
		if (flag_MSDnonAffine_t0) mkdir(MSDnonAffine_t0_dir);

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

		if (flag_CN)
		{
			config_t0.compute_CN(CN_rcut);
			config_t0.CN_to_file(CNdir + "cn.0");
		}
		vector<double> y_sum(config_t0.get_particle().size(), 0);
		vector<double> x_sum(config_t0.get_particle().size(), 0);
		vector<double> z_sum(config_t0.get_particle().size(), 0);
		for (size_t i = 0; i < config_t0.get_particle().size(); i++)
		{
			const Particle& pa = config_t0.get_particle()[i];
			x_sum[pa.id] = pa.rx;
			y_sum[pa.id] = pa.ry;
			z_sum[pa.id] = pa.rz;
		}

		size_t ncounter = 0;
		Configuration_ParticleDynamic config_last_moment;
		for (size_t imoment = start_moment; imoment <= moment_number; imoment++)
		{
			size_t istep = imoment * delta_step + t0_shift_timestep;
			ncounter++;
			string str_istep = to_string(istep);
			string config_t_fname = data_fpath + fname_prefix + str_istep + fname_postfix;
			Configuration_ParticleDynamic config_t(config_t_fname, Configuration::BoxType::tilt, shear_data_pairstyle);
			if (flag_CN)
			{
				config_t.compute_CN(CN_rcut);
				config_t.CN_to_file(CNdir + "cn." + str_istep);
			}

			if (flag_MSD)
			{
				config_t.compute_msd(config_t0);
				config_t.to_file_MSD(MSDdir + "MSD." + str_istep);
			}

			if (flag_MSDnonAffine)
			{
				config_t.compute_shear_MSDnonAffine(config_t0, shear_direction, rate);
				config_t.to_file_nonAffineMSD(MSDnonAffinedir + "MSDnonAffine." + str_istep);
			}

			if (flag_MSDnonAffine_t0)
			{
				config_t.compute_shear_MSDnonAffine_t0(config_t0, shear_direction, rate);
				config_t.to_file_nonAffineMSD(MSDnonAffine_t0_dir + "MSDnonAffine_t0." + str_istep);
			}

			if (flag_MSDnonAffine_ave_gradient)
			{
				for (size_t i = 0; i < config_t.get_particle().size(); i++)
				{
					const Particle& pa = config_t.get_particle()[i];
					x_sum[pa.id] += config_t.get_rx_real(pa);//pa.rx;
					y_sum[pa.id] += config_t.get_ry_real(pa);// pa.ry;
					z_sum[pa.id] += config_t.get_rz_real(pa);// pa.rz;
				}
				vector<double> y_ave = y_sum;
				for (size_t i = 0; i < y_ave.size(); i++)
				{
					y_ave[i] /= ncounter;
				}
				config_t.compute_shear_MSDnonAffine(config_t0, shear_direction, rate, y_ave);
				config_t.to_file_nonAffineMSD(MSDnonAffine_gradient_ave_dir + "MSDnonAffineGrdAve." + str_istep);
			}

			if (flag_computeFlowDisplacement)//!计算单粒子流场方向位移，以进一步计算流场问题，可以淘汰之前剪切带程序。
			{
				if (imoment == start_moment)
				{
					config_t.compute_shear_flow_displacement(config_t0, shear_direction);
				}
				else
				{
					config_t.compute_shear_flow_displacement(config_last_moment, shear_direction);
				}
				config_t.to_file_flow_displacement(flow_displacement_dir + "FlowDisplacement." + str_istep);
			}

			if (flag_computeFlowAveVelocity)
			{
				if (imoment == start_moment)
				{
					config_t.compute_shear_flow_ave_velocity(config_t0, shear_direction);
				}
				else
				{
					config_t.compute_shear_flow_ave_velocity(config_last_moment, shear_direction);
				}
				config_t.to_file_flow_ave_velocity(flow_ave_velocity_dir + "FlowAveVelocity." + str_istep);
			}

			bool flag_store_last_moment = flag_computeFlowAveVelocity || flag_computeFlowDisplacement;
			if (flag_store_last_moment)
			{
				config_last_moment = config_t;
				cout << config_last_moment.get_timestep() << '\n';
			}
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