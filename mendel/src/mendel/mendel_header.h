/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _MENDEL_H
#define _MENDEL_H

#include <utils/otools.h>

class mendel {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//FAM DATA
	std::vector < std::string > kids;
	std::vector < std::string > fathers;
	std::vector < std::string > mothers;

	//SAMPLE DATA
	std::vector < std::string > samples;
	std::vector < int > mothers_idx;
	std::vector < int > fathers_idx;
	std::vector < int > mendel_errors;
	std::vector < int > mendel_totals;

	//CONSTRUCTOR
	mendel();
	~mendel();

	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();
	void read_files_and_initialise();

	//
	void readPedigree(std::string fped);
	void checkMendel(int * genotypes, float & maf, int & m_errors, int & m_totals);
	void check();
	void check(std::vector < std::string > & args);
};


#endif


