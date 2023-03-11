/*******************************************************************************
 * Copyright (C) 2020 Olivier Delaneau, University of Lausanne
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

#include <mendel/mendel_header.h>

using namespace std;

void mendel::readPedigree(string fped) {
	vrb.title("Reading pedigree in [" + fped + "]");
	string buffer;
	vector < string > tokens;
	input_file fd_ped(fped);
	if (fd_ped.fail()) vrb.error("Cannot open PED file");
	while (getline(fd_ped, buffer)) {
		stb.split(buffer, tokens);
		if (tokens.size() < 3) vrb.error("Problem in pedigree file; each line should have 3 columns at least");
		kids.push_back(tokens[0]);
		fathers.push_back(tokens[1]);
		mothers.push_back(tokens[2]);
	}
	fd_ped.close();
	vrb.bullet("#families = " + stb.str(kids.size()));
}

void mendel::checkMendel(int * genotypes, float & maf, int & m_errors, int & m_totals) {

	//Get Major
	unsigned int nAC = 0, nAN = 0;
	for (int kidx = 0 ; kidx < samples.size() ; kidx++) {
		if (genotypes[2*kidx+0] != bcf_gt_missing && genotypes[2*kidx+1] != bcf_gt_missing) {
			nAC += (bcf_gt_allele(genotypes[2*kidx+0])==1) + (bcf_gt_allele(genotypes[2*kidx+1])==1);
			nAN += 2;
		}
	}
	maf = nAC * 1.0f / nAN;
	bool major = (maf > 0.5f);

	//Check Mendel
	m_errors = m_totals = 0;
	for (int kidx = 0 ; kidx < samples.size() ; kidx++) {

		int kg = -1;
		if (genotypes[2*kidx+0] != bcf_gt_missing && genotypes[2*kidx+1] != bcf_gt_missing)
			kg = (bcf_gt_allele(genotypes[2*kidx+0])==1) + (bcf_gt_allele(genotypes[2*kidx+1])==1);

		int fg = -1;
		if ((fathers_idx[kidx]>=0) && (genotypes[2*fathers_idx[kidx]+0] != bcf_gt_missing) && (genotypes[2*fathers_idx[kidx]+1] != bcf_gt_missing))
			fg = (bcf_gt_allele(genotypes[2*fathers_idx[kidx]+0])==1) + (bcf_gt_allele(genotypes[2*fathers_idx[kidx]+1])==1);

		int mg = -1;
		if ((mothers_idx[kidx]>=0) && (genotypes[2*mothers_idx[kidx]+0] != bcf_gt_missing) && (genotypes[2*mothers_idx[kidx]+1] != bcf_gt_missing))
			mg = (bcf_gt_allele(genotypes[2*mothers_idx[kidx]+0])==1) + (bcf_gt_allele(genotypes[2*mothers_idx[kidx]+1])==1) ;

		int error = 0;
		if (kg>=0 && fg>=0 && mg>=0) {
			if (fg == 0 && mg == 0 && kg == 1) { error = 1;}
			if (fg == 0 && mg == 0 && kg == 2) { error = 1;}
			if (fg == 0 && mg == 1 && kg == 2) { error = 1;}
			if (fg == 0 && mg == 2 && kg == 0) { error = 1;}
			if (fg == 0 && mg == 2 && kg == 2) { error = 1;}
			if (fg == 1 && mg == 0 && kg == 2) { error = 1;}
			if (fg == 1 && mg == 2 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 0 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 0 && kg == 2) { error = 1;}
			if (fg == 2 && mg == 1 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 2 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 2 && kg == 1) { error = 1;}
		}
		if (kg>=0 && fg>=0 && mg<0) {
			if (fg == 0 && kg == 2) { error = 1;}
			if (fg == 2 && kg == 0) { error = 1;}
		}
		if (kg>=0 && fg<0 && mg>=0) {
			if (mg == 0 && kg == 2) { error = 1;}
			if (mg == 2 && kg == 0) { error = 1;}
		}

		int total = 0;
		if (kg>=0 && fg>=0 && mg>=0) {
			//total = (!major && kg!=0 && fg!=0 && mg!=0) || (major && kg!=2 && fg!=2 && mg!=2);
			if (major) total = (kg!=2) || (fg!=2) || (mg!=2);
			else total = (kg!=0) || (fg!=0) || (mg!=0);
		}
		if (kg>=0 && fg>=0 && mg<0) {
			//total = (!major && kg!=0 && fg!=0) || (major && kg!=2 && fg!=2);
			if (major) total = (kg!=2) || (fg!=2);
			else total = (kg!=0) || (fg!=0);
		}
		if (kg>=0 && fg<0 && mg>=0) {
			//total = (!major && kg!=0 && mg!=0) || (major && kg!=2 && mg!=2);
			if (major) total = (kg!=2) || (mg!=2);
			else total = (kg!=0) || (mg!=0);
		}

		if (error > total) {
			cout << endl << "DIOGO's VERBOSE kidx=" << kidx << " af=" << maf << " nAC=" << nAC << " nAN=" << nAN << " major=" << major << " kg=" << kg << " fg=" << fg << " mg=" << mg << " error=" << error << " total=" << total << endl;
		}

		mendel_errors[kidx] += error;
		mendel_totals[kidx] += total;
		m_errors += error;
		m_totals += total;
	}
}
