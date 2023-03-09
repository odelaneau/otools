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

#include <lifter/lifter_header.h>
#include <containers/chain_file.h>

using namespace liftover;
using namespace std;

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

void lifter::lift() {
	tac.clock();
	string finput = options["input"].as < string > ();
	string foutput = options["output"].as < string > ();
	string fchain =  options["chain"].as < string > ();

	readFasta();

	vrb.title("Reading chain file in [" + fchain  + "]");
	map<std::string, Target> targets = liftover::open_chainfile(fchain );

	vrb.title("Reading data in [" + finput + "]");

	//Opening input file
	bcf_srs_t * sr =  bcf_sr_init();
	if (options["thread"].as < int > () > 1) bcf_sr_set_threads(sr, options["thread"].as < int > ());
    if (!(bcf_sr_add_reader (sr, finput.c_str()))) {
    	switch (sr->errnum) {
		case not_bgzf: vrb.error("File not compressed with bgzip!"); break;
		case idx_load_failed: vrb.error("Impossible to load index file!"); break;
		case file_type_error: vrb.error("File format not detected by htslib!"); break;
		default : vrb.error("Unknown error!");
		}
	}
    int nsamples = bcf_hdr_nsamples(sr->readers[0].header);
    vrb.bullet("#samples = " + stb.str(nsamples));

	//Opening input file
    string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (foutput.size() > 6 && foutput.substr(foutput.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (foutput.size() > 3 && foutput.substr(foutput.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(foutput.c_str(),file_format.c_str());
	if (options["thread"].as < int > () > 1) hts_set_threads(fp, options["thread"].as < int > ());
	bcf_hdr_t * hdr = sr->readers[0].header;

	if (bcf_hdr_write(fp, hdr) < 0) vrb.error("Failing to write VCF/header in [" + foutput + "]");

    // Declare arrays for data
    int ngt, ngt_arr = 0; int * gt_arr = NULL;		//Genotype

    //Read data
	bcf1_t * line_data;
	int n_parsed = 0, n_success = 0, n_nfound = 0, n_mfound = 0, n_negstrand = 0, n_refallele = 0, n_diffchr = 0;
	while(bcf_sr_next_line (sr)) {

		line_data =  bcf_sr_get_line(sr, 0);

		bcf_unpack(line_data, BCF_UN_STR);
		string chr = bcf_hdr_id2name(hdr, line_data->rid);
		int pos = line_data->pos;
		string ref = string(line_data->d.allele[0]);

		vector< Match > matches = targets[chr][pos];

		if (matches.size() == 1) {
			if (matches[0].contig == chr) {
				if (matches[0].fwd_strand) {
					int new_pos0 = matches[0].pos;
					assert(new_pos0 + ref.size() <= refseq.size());
					string new_refA = refseq.substr(new_pos0, ref.size());
					if (new_refA == ref) {
						line_data->pos = new_pos0;
						if (bcf_write1(fp, hdr, line_data) < 0) vrb.error("Failing to write VCF/record");
						n_success++;
					} else n_refallele++;
				} else n_negstrand++;
			} else n_diffchr++;
		} else if (matches.size() == 0) n_nfound++;
		else n_mfound++;
		n_parsed++;
	}
	free(gt_arr);
	bcf_destroy1(line_data);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");

	vrb.title("Writing lifted-over data in [" + foutput + "]");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF uncompressed / N=" + stb.str(nsamples) + " (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF compressed / N=" + stb.str(nsamples) + " (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF compressed / N=" + stb.str(nsamples) + " (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
	vrb.bullet("#records parsed = " + stb.str(n_parsed));
	vrb.bullet("#records successfully lifted-over = " + stb.str(n_success));
	vrb.bullet("#records NOT lifted-over = " + stb.str(n_nfound+n_mfound+n_negstrand+n_refallele+n_diffchr));
	vrb.bullet("   - position = " + stb.str(n_nfound));
	vrb.bullet("   - multi-match = " + stb.str(n_mfound));
	vrb.bullet("   - negative strand = " + stb.str(n_negstrand));
	vrb.bullet("   - unmatching REF allele = " + stb.str(n_refallele));
	vrb.bullet("   - different contig = " + stb.str(n_diffchr));

	//step2: Measure overall running time
	vrb.title("Total running time = " + stb.str(tac.abs_time()) + " seconds");

}
