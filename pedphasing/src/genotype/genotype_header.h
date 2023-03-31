#ifndef _GENOTYPE_H
#define _GENOTYPE_H

#include <utils/otools.h>

class genotype {
public:

	std::vector < std::string > vec_names;				//samples ids in std::vector
	std::map < std::string, int > map_names;				//samples ids in map
	std::vector < int > fathers;					//father ids
	std::vector < int > mothers;					//mother ids

	std::vector < std::string > chr, id, ref, alt;
	std::vector < int > pos;

	std::vector < int > mendel_error;
	std::vector < int > mendel_total;

	std::vector < std::vector < bool > > gen1;
	std::vector < std::vector < bool > > gen2;
	std::vector < std::vector < bool > > phas;
	std::vector < std::vector < bool > > miss;

	genotype();
	~genotype();

	void readPedigrees(std::string);
	void readGenotypes(std::string, std::string);
	void writeGenotypes(std::string);

	bool solveTrio(int locus, int cidx, int fidx, int midx);
	bool solveDuoFather(int locus, int cidx, int pidx);
	bool solveDuoMother(int locus, int cidx, int pidx);
	void solvePedigrees();
};

#endif
