//2020.12.16, new verison add function: custom delimiter
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

using std::string;
using std::ifstream;
using std::vector;
using std::iostream;
using std::cerr;
using std::endl;
using std::istringstream;

std::string& trim(string& s);
class Input
{
private:
	static const auto LINE_MAX = std::numeric_limits<std::streamsize>::max();  //define max line size
	ifstream infile;
	string fname; //file name
	size_t headline; //head lines
	size_t linePointer; //line position indicator
	vector<double> data; //number data
	vector<string> str_data; //string data
	size_t totalLine; //total lines of file
	size_t totalUnemptyLine;  //total umempty lines of file
public:
	Input(string _fname, size_t _headline); //constructor
	Input(); //default constructor

	inline const string& get_fname() { return fname; }; //return file name
	inline const vector<double>& get_data() const { //return number data line
		if (!data.empty())
		{
			return data;
		}
		else
		{
			cerr << "line " << linePointer << " is EMPTY" << endl;
			throw "trying get empty data";
		}

	}
	inline size_t get_linep() const { return linePointer; } //return line position
	inline size_t get_totall() const { return totalLine; } //return total line
	inline size_t get_unempty_line() const { return totalUnemptyLine; } //return total umempty line
	inline const vector<string>& get_data_str() const { return str_data; } //return str line

	inline void reset_filename(string _fname) { fname = _fname; }; //reset file name
	bool open_file(); //open file, EXCUTE THIS BEFORE READ FILE!!!
	void close_file(); //close file
	void skiphead(); //skip head lines, if you want
	size_t move_to_line(size_t _line); // move to specific line
	size_t read_line_data(char delimiter = ' ', bool skip_empty = true); // read a line of data, with custom delimiter
	size_t read_line_str(size_t _num); // read (some) line(s) of string
	size_t skip_line(size_t _num); // skip line(s)

	bool check_EOF() const; //check eof
};

