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

#include <diploidizer/diploidizer_header.h>

using namespace std;

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

void diploidizer::diploidize() {
	tac.clock();
	string finput = options["input"].as < string > ();
	string foutput = options["output"].as < string > ();
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
    int ngt_input, ngt_arr_input = 0; int * gt_arr_input = NULL;		//Input genotypes
    int * gt_arr_output = (int *)malloc(nsamples * 2 * sizeof(int));

    //Read data
	bcf1_t * line_data;
	int line_parsed = 0;
	while(bcf_sr_next_line (sr)) {

		line_data =  bcf_sr_get_line(sr, 0);

		if (line_data->n_allele == 2) {

			//read genotypes
			ngt_input = bcf_get_genotypes(hdr, line_data, &gt_arr_input, &ngt_arr_input);
			int max_ploidy = ngt_input/nsamples;
			assert(max_ploidy == 1 || max_ploidy == 2);

			for(int i = 0 ; i < nsamples ; i ++) {
				gt_arr_output[2 * i + 0] = gt_arr_input[max_ploidy * i + 0];

				if (max_ploidy == 1) {
					gt_arr_output[2 * i + 1] = gt_arr_input[max_ploidy * i + 0];
				} else if (gt_arr_input[max_ploidy * i + 1] == bcf_int32_vector_end) {
					gt_arr_output[2 * i + 1] = gt_arr_input[max_ploidy * i + 0];
				} else {
					gt_arr_output[2 * i + 1] = gt_arr_input[max_ploidy * i + 1];
				}
			}

			bcf_update_genotypes(hdr, line_data, gt_arr_output, bcf_hdr_nsamples(hdr)*2);
			if (bcf_write1(fp, hdr, line_data) < 0) vrb.error("Failing to write VCF/record");
		}
		line_parsed++;
	}
	free(gt_arr_input);
	free(gt_arr_output);
	bcf_destroy1(line_data);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(nsamples) + " / L=" + stb.str(line_parsed) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(nsamples) + " / L=" + stb.str(line_parsed) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(nsamples) + " / L=" + stb.str(line_parsed) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}
