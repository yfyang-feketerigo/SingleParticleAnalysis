// 20220825 safely cast double to integer
#include "configuration.h"
using std::cerr;
using std::endl;
Configuration::Configuration(std::string config_file, BoxType _boxtype, PairStyle _pairstyle)
{
	clog << "#LAMMPS data file reader..." << '\n';
	clog << "#Ivan Young@CIAC 20201020" << '\n';
	clog << "Reading configuration data file " << config_file << " ..." << '\n';
	//timestep = _time;
	filename = config_file;
	string firstline;
	std::ifstream infile;
	HEAD_INFO_LINE = 0;
	infile.open(config_file);
	if (infile.is_open())
	{
		clog << "Processing data file head information..." << '\n';
		//��ȡʱ�䲽��Ϣ
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

		//����������Ϣ
		infile.ignore(LINE_SKIP_MAX, '\n'); //2nd line

		infile >> particle_num;
		infile.ignore(LINE_SKIP_MAX, '\n'); //3rd line

		infile >> type_num;
		infile.ignore(LINE_SKIP_MAX, '\n'); //4th line

		infile.ignore(LINE_SKIP_MAX, '\n'); //5th line

		HEAD_INFO_LINE = 5;

		/*
		* ����x�߽�
		*/
		infile >> xlo >> xhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //6th line
		HEAD_INFO_LINE++;

		/*
		* ����y�߽�
		*/
		infile >> ylo >> yhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //7th line
		HEAD_INFO_LINE++;

		/*
		* ����z�߽�
		*/
		infile >> zlo >> zhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //8th line
		HEAD_INFO_LINE++;

		if (_boxtype == BoxType::tilt) //��б���ӵĴ���
		{
			infile >> xy >> xz >> yz;
			infile.ignore(LINE_SKIP_MAX, '\n'); //9th line
			HEAD_INFO_LINE++;
		}
		infile.ignore(LINE_SKIP_MAX, '\n'); //10th line
		HEAD_INFO_LINE++;


		string string_mass_info; //��������������Ϣ
		for (size_t i = 0; i < 3 + type_num; i++)
		{
			getline(infile, string_mass_info);
			strvec_mass_info.push_back(string_mass_info);
		} //11th - 14th line
		HEAD_INFO_LINE += (3 + type_num);
		//infile.ignore(LINE_SKIP_MAX, '\n'); //15th line

		size_t total_pair_line = 0;
		size_t pair_info_space_line = 3;

		switch (_pairstyle) //��������������Ϣ
		{
		case Configuration::PairStyle::single:
			total_pair_line = type_num;
			break;
		case Configuration::PairStyle::pair:
			total_pair_line += type_num;
			total_pair_line += type_num * (type_num - 1) / 2;
			break;
		case Configuration::PairStyle::none:
			total_pair_line = 0;
			pair_info_space_line = 0;
		default:
			break;
		}
		string string_pair_info;
		for (size_t i = 0; i < pair_info_space_line + total_pair_line; i++)
		{
			getline(infile, string_pair_info);
			strvec_pair_info.push_back(string_pair_info);
		} //15th - 19th line
		HEAD_INFO_LINE += strvec_pair_info.size();


		getline(infile, str_atoms_info); //����������Ϣ
		HEAD_INFO_LINE++;
		infile.ignore(LINE_SKIP_MAX, '\n');
		HEAD_INFO_LINE++;
		clog << "Head information has been processed" << '\n';
	}
	else
	{
		cerr << "File " << config_file << " open failed!" << endl;
		throw std::runtime_error(config_file.c_str());
	}
	infile.close();

	/*
	* ��Ļ��������data�ļ���Ϣ������У����
	*/
	clog << firstline << '\n';;
	clog << "Configuration data file " << config_file << " has " << particle_num << " particles" << '\n';
	clog << "Configuration has " << type_num << " particle type(s)" << '\n';
	clog << "Time Step: " << timestep << '\n';
	for (size_t i = 0; i < strvec_mass_info.size(); i++)
	{
		clog << strvec_mass_info[i] << '\n';
	}

	for (size_t i = 0; i < strvec_pair_info.size(); i++)
	{
		clog << strvec_pair_info[i] << '\n';
	}
	clog << '\n' << "Box Parameters:" << '\n';
	clog << "xlo, xhi: " << xlo << ' ' << xhi << '\n';
	clog << "ylo, yhi: " << ylo << ' ' << yhi << '\n';
	clog << "zlo, zhi: " << zlo << ' ' << zhi << '\n';
	clog << "xy, xz, yz: " << xy << ' ' << xz << ' ' << yz << '\n';
	//clog << "head line info read!" << endl;
	clog << '\n';
	clog << str_atoms_info << '\n';
	clog << '\n' << "File HEAD LINE: " << HEAD_INFO_LINE << '\n';
	clog << "File GAP LINE: " << GAP_LINE << '\n';
	Input::Input in_data(config_file, HEAD_INFO_LINE);
	clog << '\n';

	/*
	* ��ʼ��������������Ϣ
	*/
	clog << "Reading coordinates..." << '\n';
	in_data.open_file();
	in_data.skiphead();
	vec_particle.resize(particle_num);
	for (size_t i = 0; i < particle_num; i++)
	{
		in_data.read_line_data();
		vec_particle[i].id = std::lround(in_data.get_data()[0]);
		vec_particle[i].type = std::lround(in_data.get_data()[1]);
		vec_particle[i].rx = in_data.get_data()[2];
		vec_particle[i].ry = in_data.get_data()[3];
		vec_particle[i].rz = in_data.get_data()[4];
		vec_particle[i].box_x = std::lround(in_data.get_data()[5]);
		vec_particle[i].box_y = std::lround(in_data.get_data()[6]);
		vec_particle[i].box_z = std::lround(in_data.get_data()[7]);
	}

	in_data.skip_line(GAP_LINE); //�����������ٶȼ����
	clog << "Coordinates have been read!" << '\n';

	/*
	* ��ʼ���������ٶ���Ϣ
	*/
	clog << "Reading velocities..." << '\n';
	for (size_t i = 0; i < particle_num; i++)
	{
		in_data.read_line_data();
		Particle& p_particle = seek_id(vec_particle, (size_t)in_data.get_data()[0]);
		p_particle.vx = in_data.get_data()[1];
		p_particle.vy = in_data.get_data()[2];
		p_particle.vz = in_data.get_data()[3];
	}
	clog << "Velocities have been read!" << '\n';
	clog << "Configuration data file " << config_file << " has been read!" << '\n';
	clog << '\n';
	infile.close();
}

size_t Configuration::__add_particle(const Particle& new_pa) //�����һ������
{
	bool flag_inbox = new_pa.rx >= xlo && new_pa.rx <= xhi
		&& new_pa.ry >= ylo && new_pa.ry <= yhi
		&& new_pa.rz >= zlo && new_pa.rz <= zhi;
	if (!flag_inbox)
	{
		cerr << "new particle coordiantion is not in box!" << endl;
		throw std::runtime_error("new particle coordiantion is not in box!");
	}
	bool flag_oldtype = true;
	for (size_t i = 0; i < vec_particle.size(); i++)
	{
		flag_oldtype = flag_oldtype && (new_pa.type == vec_particle[i].type);
	}
	if (!flag_oldtype)
	{
		cerr << "WARNING: NEW particle type add, make sure pair&mass info be modified" << endl;
	}
	particle_num++;
	vec_particle.push_back(new_pa);
	return vec_particle.size();
}

void Configuration::to_data(string fname, BoxType _boxtype) //��lammps data�ļ���ʽ���
{
	ofstream ofile;
	ofile.open(fname);
	if (!ofile.is_open())
	{
		cerr << fname << " open failed" << endl;
		throw std::runtime_error((fname + " open failed").c_str());
	}
	ofile << "LAMMPS data file via C++, Configuration class, timestep = " << timestep << '\n';
	ofile << '\n';
	ofile << particle_num << " atoms" << '\n';
	ofile << type_num << " atom types" << '\n';
	ofile << '\n';

	ofile << xlo << " " << xhi << " " << "xlo " << "xhi" << '\n';
	ofile << ylo << " " << yhi << " " << "ylo " << "yhi" << '\n';
	ofile << zlo << " " << zhi << " " << "zlo " << "zhi" << '\n';
	if (_boxtype == BoxType::tilt)
	{
		ofile << xy << " " << xz << " " << yz << " " << "xy " << "xz " << "yz";
	}
	for (size_t i = 0; i < strvec_mass_info.size(); i++)
	{
		ofile << strvec_mass_info[i] << '\n';
	}
	for (size_t i = 0; i < strvec_pair_info.size(); i++)
	{
		ofile << strvec_pair_info[i] << '\n';
	}
	ofile << str_atoms_info << '\n';
	ofile << '\n';
	for (size_t i = 0; i < vec_particle.size(); i++)
	{
		Particle& pa = vec_particle[i];
		ofile << pa.id << " " << pa.type << " " << pa.rx << " " << pa.ry << " " << pa.rz << " "
			<< pa.box_x << " " << pa.box_y << " " << pa.box_z << '\n';
	}
	ofile << '\n';
	ofile << "Velocities" << '\n';
	ofile << '\n';
	for (size_t i = 0; i < vec_particle.size(); i++)
	{
		Particle& pa = vec_particle[i];
		ofile << pa.id << " " << pa.vx << " " << pa.vy << " " << pa.vz << '\n';
	}
	ofile.close();
	return;
}

void Configuration::to_dump(string ofname, string opath, string style) const //��lammps dumpo�ļ���ʽ���, style ����ָ��Ϊ��ͬ��ʽ
{
	clog << "converting data file to dump file..." << '\n';
	clog << "dump style: " << style << '\n';
	clog << "dump path: " << opath << '\n';
	clog << "dump file name: " << ofname << '\n';
	if (!boost::filesystem::exists(opath))
	{
		boost::filesystem::create_directories(opath);
	}
	string full_opath = opath + ofname;

	ofstream ofile;
	ofile.open(full_opath);
	if (!ofile.is_open())
	{
		cerr << full_opath << " open failed" << endl;
		throw std::runtime_error((full_opath + " open failed").c_str());
	}

	double lx = xhi - xlo;
	double ly = yhi - ylo;
	double lz = zhi - zlo;
	if (style == "yjruan")
	{
		ofile << "ITEM: TIMESTEP" << '\n';
		ofile << timestep << '\n';
		ofile << "ITEM: NUMBER OF ATOMS" << '\n';
		ofile << particle_num << '\n';
		ofile << "ITME: BOX BOUNDS xy xz yz pp pp pp" << '\n';
		ofile.precision(15);

		/*
		* ע��˴��߽�Ĵ���data�ļ���dump�ļ���ʽ����ͬ
		*/
		double xlo_bound = xlo + std::min({ 0.0, xy, xz, xy + xz });
		double xhi_bound = xhi + std::max({ 0.0, xy, xz, xy + xz });

		double ylo_bound = ylo + std::min({ 0.0, yz });
		double yhi_bound = yhi + std::max({ 0.0, yz });

		double zlo_bound = zlo;
		double zhi_bound = zhi;

		ofile << xlo_bound << " " << xhi_bound << " " << xy << '\n';
		ofile << ylo_bound << " " << yhi_bound << " " << xz << '\n';
		ofile << zlo_bound << " " << zhi_bound << " " << yz << '\n';
		ofile << "ITEM: ATOMES id xu yu zu" << '\n';
		for (size_t i = 0; i < vec_particle.size(); i++)
		{
			const Particle& p_pa = vec_particle[i];
			double xu = get_rx_real(p_pa);
			double yu = get_ry_real(p_pa);
			double zu = get_rz_real(p_pa);
			size_t id = p_pa.id;
			ofile << id << " ";
			ofile << xu << " ";
			ofile << yu << " ";
			ofile << zu << " ";
			ofile << '\n';
		}
		return;
	}

	cerr << "wrong dump style: " << style << '\n';
	throw std::runtime_error(("wrong dump style: " + style).c_str());

}