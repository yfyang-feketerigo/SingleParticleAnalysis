//2020.12.16, new verison add function: custom delimiter
//2021.11.29, 规范namespace
//2022.3.22, 移除自动统计文件行数功能
/*
* 用于按行读取数据块
* 每次读取一行数据，以vecotr<double>形式存储，返回
* 每行数据长度可以不一致
* 可以按行读取字符
* 注意分割符
*/
#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <limits>
#include <utility>
#include <exception>
#include <stdexcept>
//#define NO_THROW_EMPTY_LINE
//#define COUNT_FILE_TOTAL_LINE
namespace Input {
	std::string& trim(std::string& s);
	std::string add_quote(const std::string& str);
	class Input
	{
	private:
		static const auto LINE_MAX = std::numeric_limits<std::streamsize>::max();  //define max line size
		std::ifstream infile;
		std::string fname; //file name
		size_t headline; //head lines
		size_t linePointer; //line position indicator
		std::vector<double> data; //number data
		std::vector<std::string> str_data; //string data
#ifdef COUNT_FILE_TOTAL_LINE
		size_t totalLine; //total lines of file
		size_t totalUnemptyLine;  //total unempty lines of file
#endif // COUNT_FILE_TOTAL_LINE
	public:
		Input(std::string _fname, size_t _headline); //constructor
		Input(); //default constructor
		Input(const Input&) = delete;
		Input& operator=(const Input&) = delete;
		inline const std::string& get_fname() { return fname; }; //return file name
		inline const std::vector<double>& get_data() const { //return number data line
			if (!data.empty())
			{
				return data;
			}
			else
			{
				std::cerr << "line " << linePointer << " is EMPTY" << std::endl;
				throw std::runtime_error("trying get empty data");
			}

		}
		inline size_t get_linep() const { return linePointer; } //return line position
#ifdef COUNT_FILE_TOTAL_LINE
		inline size_t get_totall() const { return totalLine; } //return total line
		inline size_t get_unempty_line() const { return totalUnemptyLine; } //return total umempty line
#endif // COUNT_FILE_TOTAL_LINE


		inline const std::vector<std::string>& get_data_str() const { return str_data; } //return str line

		inline void reset_filename(std::string _fname) { fname = _fname; }; //reset file name
		bool open_file(); //open file, EXCUTE THIS BEFORE READ FILE!!!
		void close_file(); //close file
		void skiphead(); //skip head lines, you will need this in most case :)
		size_t move_to_line(size_t _line); // move to specific line
		size_t read_line_data(char delimiter = ' ', bool skip_empty = true); // read a line of data, with custom delimiter
		size_t read_line_str(size_t _num = 1); // read (some) line(s) of string
		size_t skip_line(size_t _num = 1); // skip line(s)

		bool check_EOF() const; //check eof
	};
}
