// Extracts random patterns from a file

#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>

// sdsl bit vector functions
#include <sdsl/bit_vectors.hpp>
// memory counter
#include "external/malloc_count/malloc_count.h"

int main (int argc, char **argv)
{

	FILE *ofile;

	if (argc < 6)
	{
		std::cerr <<
			 "Usage: genpatterns <file> <length> <number> <patterns file> <circular> <verbose>\n" <<
			 "  randomly extracts <number> substrings of length <length> from <file>,\n" <<
			 "  The output file, <patterns file> is a fasta file with one pattern per line\n" <<
			 "  <circular> is a flag, if circular=0 it gives in output only linear pattern,\n  if circular=1 both linear and circular patterns are computed.\n" <<
			 "  <verbose> is a flag, if verbose=0 no message is written,\n if verbose=1 the verbose mode is switched one\n";
		exit(1);
	}

	srand(261222); 

	// read input parameters
	std::string ifile_name = argv[1];
	int64_t plen = atoi(argv[2]);
	int64_t npat = atoi(argv[3]);
	std::string ofile_name = argv[4];
	bool circular = atoi(argv[5]);
	bool verbose = atoi(argv[6]);
	int64_t no_circ = 0;

	// print input parameters
	if( verbose )
	{
		{
			std::cout << "Input file: "    << ifile_name << std::endl; 
			std::cout << "Patterns length: "    << plen  << std::endl; 
			std::cout << "No. patterns: "       << npat  << std::endl; 
			std::cout << "output file: "  << ofile_name  << std::endl;
		}
		if( circular )
		{
			std::cout << "Computing both linear and circular patterns"  << std::endl;
		}
		else
		{
			std::cout << "Computing only linear patterns"  << std::endl;
		}
	}

	// compute string bounduaries bit vectors
	std::vector<int64_t> onset_beg;
	std::vector<int64_t> onset_end;
	std::ifstream input(ifile_name);
	std::string line, DNA_sequence;
    int64_t sum = 0, ns = 0;
    while(std::getline(input, line)) {
        // skip empty lines
        if(line.empty()){ continue; }
        // header of a new sequence
        if (line[0] == '>') {
        	ns++;
            size_t seqlen = DNA_sequence.size();
            if(seqlen > 0 ) sum += seqlen + 1;
            if(seqlen > 0){ onset_end.push_back(sum-2); }
            sum += line.size()+1;
            onset_beg.push_back(sum);
            // insert previous DNA sequence
            // the current DNA sequence
            DNA_sequence.clear();
        }
        else {
            // add new line
            DNA_sequence += line;
        }
    }
    // add last sequence
    sum += DNA_sequence.size() + 1;
    onset_end.push_back(sum-2);
    input.clear();

    // init buffer
    char * buffer = new char [plen];
    // compute select data structures
    // beginning positions
    sdsl::sd_vector_builder builder_beg(sum,onset_beg.size());
    for(auto idx: onset_beg){ builder_beg.set(idx); }
    sdsl::sd_vector<> bv_beg = sdsl::sd_vector<>(builder_beg);
	sdsl::sd_vector<>::select_1_type s_beg = sdsl::sd_vector<>::select_1_type(&bv_beg);
	onset_beg.clear();
	// ending positions
    sdsl::sd_vector_builder builder_end(sum,onset_end.size());
    for(auto idx: onset_end){ builder_end.set(idx); }
    sdsl::sd_vector<> bv_end = sdsl::sd_vector<>(builder_end);
	sdsl::sd_vector<>::select_1_type s_end = sdsl::sd_vector<>::select_1_type(&bv_end);
	onset_end.clear();

	// initialize output file
	ofile = fopen(ofile_name.c_str(), "w+");
	int64_t pat_id = 0;
	while( npat > 0 )
	{
		// select sequence
		int64_t nseq = rand()%ns + 1;
		// take sequence length
		int64_t beginning = s_beg(nseq);
		int64_t ending =    s_end(nseq);
		// sampling position
		int64_t pat_pos = beginning + rand()%(ending - beginning);
		// check pattern circularity
		if( (ending - pat_pos) < plen )
		{
			if( !circular || (ending - beginning + 1) < plen ){ continue; }
			input.seekg(pat_pos, std::ios::beg);
			input.read(buffer,(ending - pat_pos + 1));
			input.seekg(beginning, std::ios::beg);
			input.read(&buffer[(ending - pat_pos + 1)],plen-(ending - pat_pos + 1));
			no_circ++;
		}
		else
		{	
			input.seekg(pat_pos, std::ios::beg);
			input.read(buffer,plen);
		}
		pat_id++;
		std::string header = (std::string)">Pattern" + std::to_string(pat_id) + (std::string)"\n";
		fwrite (&header[0] , sizeof(char), header.size(), ofile);
		//ofile << header;
		fwrite (buffer , sizeof(char), plen, ofile);
		putc('\n', ofile);
		npat--;
	}

	std::cout << "Done! number of sampled circular pattern: " << no_circ << "\n";

	input.close();
	fclose(ofile);

	return 0;
}
