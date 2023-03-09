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

void mendel::check() {
	tac.clock();

	//
	readPedigree(options["pedigree"].as < string > ());

	//
	string finput = options["input"].as < string > ();
	string foutput = options["output"].as < string > ();
	string region = options["region"].as < string > ();
	vrb.title("Reading data in [" + finput + "]");

	//Opening input file
	bcf_srs_t * sr =  bcf_sr_init();
	if (options["thread"].as < int > () > 1) bcf_sr_set_threads(sr, options["thread"].as < int > ());
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "]");
	if (!(bcf_sr_add_reader (sr, finput.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf: vrb.error("File not compressed with bgzip!"); break;
		case idx_load_failed: vrb.error("Impossible to load index file!"); break;
		case file_type_error: vrb.error("File format not detected by htslib!"); break;
		default : vrb.error("Unknown error!");
		}
	}

    //Sample IDs
    int nsamples = bcf_hdr_nsamples(sr->readers[0].header);
    vrb.bullet("#samples = " + stb.str(nsamples));
    fathers_idx = vector < int > (nsamples, -1);
    mothers_idx = vector < int > (nsamples, -1);
    mendel_errors = vector < int > (nsamples, 0);
	mendel_totals = vector < int > (nsamples, 0);

    map < std::string , int > mapF2S;
    for (int i = 0 ; i < nsamples ; i ++) {
    	samples.push_back(std::string(sr->readers[0].header->samples[i]));
    	mapF2S.insert(pair < std::string, int > (samples.back(), i));
    }
    unsigned int ntrios = 0, nduosF = 0, nduosM = 0;
    for (int f = 0 ; f < kids.size() ; f ++) {
    	map < std::string , int > ::iterator itK = mapF2S.find(kids[f]);
    	map < std::string , int > ::iterator itF = mapF2S.find(fathers[f]);
    	map < std::string , int > ::iterator itM = mapF2S.find(mothers[f]);

    	if (itK != mapF2S.end() && itF != mapF2S.end() && itM != mapF2S.end()) {
    		fathers_idx[itK->second] = itF->second;
    		mothers_idx[itK->second] = itM->second;
    		ntrios++;
    	}

    	if (itK != mapF2S.end() && itF != mapF2S.end() && itM == mapF2S.end()) {
    		fathers_idx[itK->second] = itF->second;
    		nduosF++;
    	}

    	if (itK != mapF2S.end() && itF == mapF2S.end() && itM != mapF2S.end()) {
    		mothers_idx[itK->second] = itM->second;
    		nduosM++;
    	}
    }
    vrb.bullet("#trios = " + stb.str(ntrios) + " | #duos_paternal = " + stb.str(nduosF) + " | #duos_maternal = " + stb.str(nduosM));

    //Read data and output to file
    output_file fdv(foutput + ".var.txt.gz");
    int ngt, ngt_arr = 0; int * gt_arr = NULL, line = 0;
    bcf1_t * line_data;
	while(bcf_sr_next_line (sr)) {
		line_data =  bcf_sr_get_line(sr, 0);
		if (line_data && line_data->n_allele == 2) {
			std::string chr = bcf_hdr_id2name(sr->readers[0].header, line_data->rid);
			int pos = line_data->pos + 1;
			std::string id = std::string(line_data->d.id);
			std::string ref = std::string(line_data->d.allele[0]);
			std::string alt = std::string(line_data->d.allele[1]);
			ngt = bcf_get_genotypes(sr->readers[0].header, line_data, &gt_arr, &ngt_arr);
			assert(ngt == 2 * nsamples);
			float maf;
			int v_errors = 0, v_totals = 0;
			checkMendel(gt_arr, maf, v_errors, v_totals);
			fdv << chr << "\t" << pos << "\t" << ref << "\t" << alt << "\t" << maf << "\t" << v_errors << "\t" << v_totals << "\t" << stb.str(v_totals?(v_errors*100.0f/v_totals):0, 2) << endl;
		}
		line++;
		if (line % 10000 == 0) vrb.bullet("Processing VCF record: [" + stb.str(line) + "]");
	}
	fdv.close();
	free(gt_arr);
	bcf_sr_destroy(sr);

    //Per sample summary
	vrb.title("Writing per sample summary in [" + foutput + "]");
	output_file fds(foutput + ".ind.txt.gz");
	for (int kidx = 0 ; kidx < samples.size() ; kidx++) {
		fds << samples[kidx] << "\t" << ((fathers_idx[kidx]>=0)?samples[fathers_idx[kidx]]:"NA") << "\t" << ((mothers_idx[kidx]>=0)?samples[mothers_idx[kidx]]:"NA") << "\t" << mendel_errors[kidx] << "\t" << mendel_totals[kidx] << endl;
	}
	fds.close();

	//step2: Measure overall running time
	vrb.title("Total running time = " + stb.str(tac.abs_time()) + " seconds");

}
