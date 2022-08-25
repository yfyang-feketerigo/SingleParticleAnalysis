//2020.12.16, new verison add function: custom delimiter
//2021.11.29, 规范namespace
//2022.2.24 规范异常抛出
#include "input.h"
namespace Input
{
	using std::string;
	using std::ifstream;
	using std::vector;
	using std::iostream;
	using std::cerr;
	using std::endl;
	using std::istringstream;

	Input::Input(string _fname, size_t _headline)
	{
		fname = _fname;
		headline = _headline;
		linePointer = 0;
#ifdef COUNT_FILE_TOTAL_LINE
		infile.open(fname);
		totalLine = 0;
		string line;
		if (infile.is_open())
		{

			while (!infile.eof())
			{

				// 打开后统计文件总行数

				totalLine += 1;
				getline(infile, line);
				if (!line.empty())
				{
					totalUnemptyLine += 1;
				}
			}

			infile.close();
		}
		else
		{
			string err_msg = "Input(): open input file " + add_quote(fname) + " failed!";
			throw std::runtime_error(err_msg.c_str());
		}
#endif // COUNT_FILE_TOTAL_LINE
	}

	Input::Input()
	{
		/*
		* 默认构造函数
		*/
		fname = "default";
		headline = 0;
		linePointer = 0;
		infile.open(fname);
#ifdef COUNT_FILE_TOTAL_LINE
		totalLine = 0;
		totalUnemptyLine = 0;
#endif // COUNT_FILE_TOTAL_LINE
	}

	bool Input::open_file()
	{
		/*
		* 打开文件
		*/
		if (!infile.is_open())
		{
			infile.open(fname);
		}
		else
		{
			//cerr << "input file has already been opened!" << endl;
			string err_msg = "open_file(): input file " + add_quote(fname) + " has already been opened!";
			throw std::runtime_error(err_msg.c_str());
		}
		return infile.is_open();
	}

	void Input::close_file()
	{
		/*
		* 关闭文件
		*/
		if (infile.is_open())
		{
			infile.close();
			linePointer = 0;
#ifdef COUNT_FILE_TOTAL_LINE
			totalLine = 0;
#endif // COUNT_FILE_TOTAL_LINE
			data.clear();
		}
		else
		{
			//cerr << "input file is not opened!" << endl;
			string err_msg = "close_file(): nput file " + add_quote(fname) + " is not opened!";
			throw std::runtime_error(err_msg.c_str());
		}
	}

	void Input::skiphead()
	{
		/*
		* 跳过文件头
		*/
		if (infile.is_open())
		{
			if (linePointer == 0)
			{
				for (size_t i = 0; i < headline; i++)
				{
					infile.ignore(LINE_MAX, '\n');
					linePointer += 1;
				}
			}
			else
			{
				//cerr << "head lines in file has already been skipped!" << endl;
				string err_msg = "skiphead(): head lines in file " + add_quote(fname) + " has already been skipped!";
				throw std::runtime_error(err_msg.c_str());
			}
		}
		else
		{
			//cerr << fname << " not open!" << endl;
			throw std::runtime_error(("skiphead(): " + add_quote(fname) + " not open!").c_str());
		}
	}

	size_t Input::move_to_line(size_t _line)
	{
		/*
		* 移动到某行
		*/
#ifdef COUNT_FILE_TOTAL_LINE
		if (_line > totalLine)
		{
			//cerr << "out of file line range!" << endl;
			//return linePointer;
			string err_msg = "move_to_line(): out of file line range: " + add_quote(fname);
			throw std::runtime_error(err_msg.c_str());
		}
#endif // COUNT_FILE_TOTAL_LINE
		infile.close();
		infile.open(fname);
		for (size_t i = 1; i < _line; i++)
		{
			infile.ignore(LINE_MAX, '\n');
		}
		return linePointer = _line;
	}

	size_t Input::read_line_data(char delimiter, bool skip_empty)
	{
		/*
		* 读取一行数据，注意是数字数据
		*/
		data.clear();
		string line;
		if (infile.is_open())
		{
			if (infile.eof())
			{
				string err_msg = "read_line_data(): file " + add_quote(fname) + " reach eof! data remain unchanged!";
				throw std::runtime_error(err_msg.c_str());
			}
			else
			{
				do
				{
					getline(infile, line);
					if (skip_empty)
					{
						if (line.empty())
						{
#ifdef NO_THROW_EMPTY_LINE
							cerr << "!!!WARNING!!!: skip reading EMPTY data line " << linePointer << " at" << fname << '\n'
								<< "this may cause SEVERE problems when using date\n";
#else
							throw std::runtime_error(("read_line_data(): trying reading EMPTY data line: " + std::to_string(linePointer) + " at file: " + add_quote(fname)).c_str());
#endif // NO_THROW_EMPTY_LINE
							linePointer += 1;
						}
					}
				} while (line.empty() && !infile.eof());
				istringstream ss(line);
				string str_num;
				while (std::getline(ss, str_num, delimiter))
				{
					str_num = trim(str_num);
					data.push_back(std::stod(str_num));
				}

			}
			linePointer += 1;
			return data.size();
		}

		else
		{
			//cerr << "input file is not open!" << endl;
			string err_msg = "read_line_data(): input file " + add_quote(fname) + " is not open!";
			throw std::runtime_error(err_msg.c_str());
		}
	}

	size_t Input::read_line_str(size_t _num)
	{
		/*
		* 读取若干行字符串
		*/
		string str_temp;
		if (!str_data.empty())
		{
			str_data.clear();
		}
		if (infile.is_open())
		{
			for (size_t i = 0; i < _num; i++)
			{
				getline(infile, str_temp);
				str_data.push_back(str_temp);
				linePointer++;
			}
		}
		return str_data.size();
	}

	size_t Input::skip_line(size_t _num)
	{
		/*
		* 跳过若干行
		*/
		if (infile.is_open())
		{
			for (size_t i = 0; i < _num; i++)
			{
				infile.ignore(LINE_MAX, '\n');
				linePointer++;
			}
		}
		else
		{
			string err_msg = "skip_lin(): input file " + add_quote(fname) + " is not open!";
			throw std::runtime_error(err_msg.c_str());
		}
		return linePointer;
	}

	bool Input::check_EOF() const
	{
		if (infile.eof())
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	string& trim(string& s)
	{
		/*
		* 裁剪字符串，去除首尾空格
		*/
		if (s.empty())
			return s;
		s.erase(0, s.find_first_not_of(' '));
		s.erase(s.find_last_not_of(' ') + 1);
		return s;
	}
	std::string add_quote(const std::string& str)
	{
		auto quoted_string = "\"" + str + "\"";
		return quoted_string;
	}
}