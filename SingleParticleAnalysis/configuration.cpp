#include "configuration.h"

size_t Configuration::GAP_LINE = 3;
Configuration::Configuration(std::string config_file, BoxType _boxtype, PairStyle _pairstyle)
{
	clog << "#LAMMPS data file reader..." << endl;
	clog << "#Ivan Young@CIAC 20201020" << endl;
	clog << "Reading configuration data file " << config_file << " ..." << endl;
	//timestep = _time;
	filename = config_file;
	string firstline;
	ifstream infile;
	HEAD_INFO_LINE = 0;
	infile.open(config_file);
	if (infile.is_open())
	{
		clog << "Processing data file head information..." << endl;
		getline(infile, firstline); //1st line
		string str_timestep;
		for (size_t i = firstline.size() - 1; i > 0; i--)
		{
			if (' ' != firstline[i])
			{
				str_timestep.insert(str_timestep.begin(), firstline[i]);
			}
			else
			{
				break;
			}
		}
		timestep = stoull(str_timestep); //string to unsigned long long
		//cout << timestep << endl;

		infile.ignore(LINE_SKIP_MAX, '\n'); //2nd line

		infile >> particle_num;
		infile.ignore(LINE_SKIP_MAX, '\n'); //3rd line

		infile >> type_num;
		infile.ignore(LINE_SKIP_MAX, '\n'); //4th line

		infile.ignore(LINE_SKIP_MAX, '\n'); //5th line

		HEAD_INFO_LINE = 5;

		infile >> xlo >> xhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //6th line
		HEAD_INFO_LINE++;

		infile >> ylo >> yhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //7th line
		HEAD_INFO_LINE++;

		infile >> zlo >> zhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //8th line
		HEAD_INFO_LINE++;

		if (_boxtype == BoxType::tilt)
		{
			infile >> xy >> xz >> yz;
			infile.ignore(LINE_SKIP_MAX, '\n'); //9th line
			HEAD_INFO_LINE++;
		}
		infile.ignore(LINE_SKIP_MAX, '\n'); //10th line
		HEAD_INFO_LINE++;


		string string_mass_info;
		for (size_t i = 0; i < 3 + type_num; i++)
		{
			getline(infile, string_mass_info);
			strvec_mass_info.push_back(string_mass_info);
		} //11th - 14th line
		HEAD_INFO_LINE += (3 + type_num);
		//infile.ignore(LINE_SKIP_MAX, '\n'); //15th line

		size_t total_pair_line = 0;
		switch (_pairstyle)
		{
		case Configuration::PairStyle::single:
			total_pair_line = type_num;
			break;
		case Configuration::PairStyle::pair:
			total_pair_line += type_num;
			total_pair_line += type_num * (type_num - 1) / 2;
			break;
		default:
			break;
		}
		string string_pair_info;
		for (size_t i = 0; i < 3 + total_pair_line; i++)
		{
			getline(infile, string_pair_info);
			strvec_pair_info.push_back(string_pair_info);
		} //15th - 19th line
		HEAD_INFO_LINE += strvec_pair_info.size();


		getline(infile, str_atoms_info);
		HEAD_INFO_LINE++;
		infile.ignore(LINE_SKIP_MAX, '\n');
		HEAD_INFO_LINE++;
		clog << "Head information has been processed" << endl;
	}
	else
	{
		cerr << "File " << config_file << " open failed!" << endl;
		throw config_file;
	}
	infile.close();

	clog << firstline << endl;;
	clog << "Configuration data file " << config_file << " has " << particle_num << " particles" << endl;
	clog << "Configuration has " << type_num << " particle type(s)" << endl;
	clog << "Time Step: " << timestep << endl;;
	for (size_t i = 0; i < strvec_mass_info.size(); i++)
	{
		clog << strvec_mass_info[i] << endl;
	}

	for (size_t i = 0; i < strvec_pair_info.size(); i++)
	{
		clog << strvec_pair_info[i] << endl;
	}
	clog << endl << "Box Parameters:" << endl;
	clog << "xlo, xhi: " << xlo << ' ' << xhi << endl;
	clog << "ylo, yhi: " << ylo << ' ' << yhi << endl;
	clog << "zlo, zhi: " << zlo << ' ' << zhi << endl;
	clog << "xy, xz, yz: " << xy << ' ' << xz << ' ' << yz << endl;
	//clog << "head line info read!" << endl;
	clog << endl;
	clog << str_atoms_info << endl;
	clog << endl << "File HEAD LINE: " << HEAD_INFO_LINE << endl;
	clog << "File GAP LINE: " << GAP_LINE << endl;
	Input in_data(config_file, HEAD_INFO_LINE);
	clog << endl;
	clog << "Reading coordinates..." << endl;
	in_data.open_file();
	in_data.skiphead();
	vec_particle.resize(particle_num);
	for (size_t i = 0; i < particle_num; i++)
	{
		in_data.read_line_data();
		vec_particle[i].id = (size_t)in_data.get_data()[0];
		vec_particle[i].type = (int)in_data.get_data()[1];
		vec_particle[i].rx = in_data.get_data()[2];
		vec_particle[i].ry = in_data.get_data()[3];
		vec_particle[i].rz = in_data.get_data()[4];
		vec_particle[i].box_x = (int)in_data.get_data()[5];
		vec_particle[i].box_y = (int)in_data.get_data()[6];
		vec_particle[i].box_z = (int)in_data.get_data()[7];
	}

	in_data.skip_line(GAP_LINE);
	clog << "Coordinates have been read!" << endl;
	clog << "Reading velocities..." << endl;
	Particle* p_particle = nullptr;
	for (size_t i = 0; i < particle_num; i++)
	{
		in_data.read_line_data();
		p_particle = seek_id(vec_particle, (size_t)in_data.get_data()[0]);
		(*p_particle).vx = in_data.get_data()[1];
		(*p_particle).vy = in_data.get_data()[2];
		(*p_particle).vz = in_data.get_data()[3];
	}
	p_particle = nullptr;
	clog << "Velocities have been read!" << endl;
	clog << "Configuration data file " << config_file << " has been read!" << endl;
	clog << endl;
	infile.close();
}
