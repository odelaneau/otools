////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <lifter/lifter_header.h>

using namespace std;

lifter::lifter() {
}

lifter::~lifter() {
}

void lifter::lift(vector < string > & args) {
	declare_options();
	parse_command_line(args);
	check_options();
	verbose_files();
	verbose_options();
	lift();
}

void lifter::readFasta() {
	tac.clock();
	string buffer;
	vector < string > tokens;
	string ffasta =  options["fasta"].as < string > ();
	string schrom =  options["chr"].as < string > ();
	bool chr_found = false;

	vrb.title("Reading fasta file in [" + ffasta  + "]");
	input_file fd(ffasta);

	while (getline(fd, buffer)) {
		if (buffer[0] == '>' && chr_found) break;
		else if (buffer[0] == '>') {
			stb.split(buffer, tokens);
			string contig = tokens[0].substr(1);
			if (contig == schrom) chr_found = true;
		} else if (chr_found) refseq+=buffer;
	}
	fd.close();
	vrb.bullet("L=" + stb.str(refseq.size()));
	vrb.bullet("Time = " + stb.str(tac.rel_time()*0.001, 2) + "s");
}
