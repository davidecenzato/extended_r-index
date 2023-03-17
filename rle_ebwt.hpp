/*
 * Implementation of RLFM-index for the eBWT.
 * 
 * This code is adapted from https://github.com/nicolaprezza/r-index.git
 *
 */

#ifndef RLE_EBWT_HPP_
#define RLE_EBWT_HPP_

#include <sdsl/wavelet_trees.hpp>
#include "sd_vector.hpp"

class rle_ebwt{

public:
	// wavelet tree
	sdsl::wt_huff<> bwt_heads;
	// BWT C vector (F column)
	std::vector<uint_t> C;
	// present characters
	//cstd::vector<bool> C_p;
	// empty constructor
	rle_ebwt(){}
	/*
	* Constructor that takes in input the files with the rle eBWT and the block size
	* and construct the data structure for rank and select queries on the rle eBWT
	*/
	rle_ebwt(std::string &headfile, std::string &lenfile, uint_t B_, bool verbose = false){
		// read heads file
		read_file(headfile.c_str(),heads);
		// read head lengths file
		read_file(lenfile.c_str(),lens);

		assert(heads.size() == lens.size());
		B = B_;

		// initialize bitvector data structures
		letter_bv.resize(128);
		// initialize onset vectors
		std::vector<uint_t> onset_main;
		std::vector< std::vector<uint_t> > onset_letter;
		onset_letter.resize(128);
		// initialize C vector and R
		C.resize(128);
		R=0;
		// iterate over heads
		for(size_t i=0;i<heads.size();++i){
			// if the run contains at lest two characters
			if(lens[i] > 1){
				// increase C vector entry
				C[heads[i]] += lens[i]-1;
				// insert lens[i]-1 0s
				BWTlength += lens[i]-1;
			}
			//runs_per_letter_bv[heads[i]].push_back(true);
			onset_letter[heads[i]].push_back(C[heads[i]]);
			//push back a bit set only at the end of a block
			if(R%B==B-1){ onset_main.push_back(BWTlength); }
			//runs_bv.push_back(R%B==B-1);
			// increase BWT length
			BWTlength++;
			// increase char counter
			C[heads[i]]++;
			// increase R
			R++;
		}
		// print stats
		if(verbose)
		{
			std::cout << "eBWT length: " << BWTlength;
			std::cout << "\nNumber of eBWT equal-letter sampled runs: " << R << std::endl;
		}
	    // check BWT length
	    #if M64
	        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
	        if(BWTlength > pow(2,32) - 1){ 
	            // the input file is too big
	            std::cerr << "Error, the eBWT length is > 4.29 GB, please use ./er-index64. exiting..." << std::endl;
	            exit(-1);
	        }
	    #endif 
		// construct the main compressed bitvector
		main_bv = sd_vector(onset_main, BWTlength);
		// construct the compressed bitvector for each character
		for(int i=0; i<128; ++i){
			// if we have at least one character
			if(C[i] > 0){
				// construct compressed bit vector
				letter_bv[i] = sd_vector(onset_letter[i], C[i]);
			}
		}
		// construct C vector
		memmove(&C[1], &C[0], 127*sizeof(uint_t));
		C[0] = 0;
		for(int i=1; i<128; ++i){ C[i] += C[i-1]; }
		// free memory
		heads.clear();
		lens.clear();
		// construct the wavalet tree for the eBWT heads
		sdsl::construct(bwt_heads, headfile.c_str(), 1);
	}

	// 2nd constructor
	rle_ebwt(std::ifstream& headfile, std::ifstream& lenfile, std::string &headstr, uint_t B_, int isize, bool verbose = false){
		// set block size
		B = B_;
		// get no runs
		headfile.seekg(0, std::ios::end);
		R = headfile.tellg();
    	headfile.seekg(0, std::ios::beg);
		// initialize bitvector data structures
		letter_bv.resize(128);
		// initialize main onset vector
		std::vector<uint_t> onset_main;
		onset_main.reserve(R/B);
		// initialize letters onset vector
		std::vector< std::vector<uint_t> > onset_letter;
		onset_letter.resize(128);
		// initialize C vector and R
		C.resize(128);
		// iterate over heads
		uint64_t currLen = 0, BWTlen = 0;
		char currHead = 0;
		// iterate over head vector
		for(size_t i=0;i<R;++i){
			// get current run length
			lenfile.read(reinterpret_cast<char*>(&currLen), isize);
			headfile >> currHead;
			// if the run contains at least two characters
			if(currLen > 1){
				// increase C vector entry
				C[currHead] += currLen-1;
				// insert lens[i]-1 0s
				BWTlength += currLen-1;
				BWTlen += currLen-1;
			}
			// create onset vector for letter bitvector
			onset_letter[currHead].push_back(C[currHead]);
			//push back a bit set only at the end of a block
			if(i%B==B-1){ onset_main.push_back(BWTlength); }
			// increase BWT length
			BWTlength++;
			BWTlen++;
			// increase char counter
			C[currHead]++;
		}
		// print stats
		if(verbose)
		{
			std::cout << "eBWT length: " << BWTlen;
			std::cout << "\nNumber of eBWT equal-letter sampled runs: " << R << std::endl;
		}
	    // check BWT length
	    #if M64 == 0
	    	// std::cout << "entra";
	        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
	        if(BWTlen > pow(2,32) - 1){ 
	            // the input file is too big
	            std::cerr << "Error, the eBWT length is > 4.29 GB, please use ./er-index64. exiting..." << std::endl;
	            exit(-1);
	        }
	    #endif 
		// construct the main compressed bitvector
		main_bv = sd_vector(onset_main, BWTlength);
		// construct the compressed bitvector for each character
		for(int i=0; i<128; ++i){
			// if we have at least one character
			if(C[i] > 0){
				// construct compressed bit vector
				letter_bv[i] = sd_vector(onset_letter[i], C[i]);
			}
		}
		// construct C vector
		memmove(&C[1], &C[0], 127*sizeof(uint_t));
		C[0] = 0;
		for(int i=1; i<128; ++i){ C[i] += C[i-1]; }
		// close streams
		headfile.close();
		lenfile.close();
		// construct the wavalet tree for the eBWT heads
		sdsl::construct(bwt_heads, headstr.c_str(), 1);
	}

	/*
	* return eBWT size
	*/
	uint_t size(){
		return BWTlength;
	}

	/*
	* return eBWT number of runs
	*/
	uint_t nrun(){
		return R;
	}

	/*
	* return a eBWT position
	*/
	char operator[](uint_t i){
		return bwt_heads[run_of(i).first];
	}

	/*
	 * return in which run the index i is contained
	 */
	uint_t run_of_position(uint_t i){

		uint_t last_block = main_bv.rank1(i);
		uint_t current_run = last_block*B;

		//current position in the string: the first of a block
		uint_t pos = 0;
		if(last_block>0){ pos = main_bv.select1(last_block-1)+1; }

		assert(pos <= i);

		//otherwise, scan at most B runs
		while(pos < i)
		{
			pos += run_at(current_run);
			current_run++;
		}

		if(pos>i) current_run--;

		//position i is inside run current_run
		assert(current_run<R);

		return current_run;

	}

	/*
	 *  <j=run of position i, last position of j-th run>
	 */
	std::pair<uint_t,uint_t> run_of(uint_t i){

		// compute current run position
		uint_t last_block = main_bv.rank1(i);
		uint_t current_run = last_block*B;

		//current position in the string: the first of a block
		uint_t pos = 0;
		if(last_block>0){ pos = main_bv.select1(last_block-1)+1; }

		// while the position is smaller than i
		while(pos < i){
			// add run length
 			pos += run_at(current_run);
			current_run++;

		}
		// if pos == i add another run length
		if(pos>i){
			current_run--;
		}else{
			pos += run_at(current_run);
		}
		return {current_run,pos-1};

	}

	/*
	 * length of i-th run
	 */
	uint_t run_at(uint_t i){
		// current head
		char c = bwt_heads[i];
		// return gap
		return letter_bv[c].gapAt(bwt_heads.rank(i,c));
	}

	/*
	 * construct rank select data structures for all bitvectors
	 */
	void construct_rank_select_dt(){
		// rank select for main bitvector
		main_bv.construct_rank_ds();
		main_bv.construct_select_ds();
		// rank select for letter bitvectors
		for(int i=0; i<128; ++i){
			if(letter_bv[i].size() > 0){
				letter_bv[i].construct_rank_ds();
				letter_bv[i].construct_select_ds();
			}
		}
	}

	/*
	 * number of c before position i
	 */
	uint_t rank(uint_t i, char c, uint_t B){
		// if c is not in the text
		// if(letter_bv[c].size()==0) return 0;
		// if i is equal the size of the eBWT
		if(i==BWTlength) return letter_bv[c].size();
		// get current run
		uint_t last_block = main_bv.rank1(i);
		uint_t current_run = last_block*B;
		// get first position of the previous block
		uint_t pos = 0;
		if( last_block>0 ){ pos = main_bv.select1(last_block-1)+1; }
		// get distance between i and previous block
		// assert(pos <= i);
		uint_t dist = i-pos;
		//otherwise, scan at most B runs
		while(pos < i){
			// get current run length
			//pos += run_length(current_run);
			pos += run_at(current_run);
			current_run++;
			// update the distance until we get to the
			// correct run
			if(pos<=i) dist = i-pos;
		}
		// get the correct run counter
		if(pos>i) current_run--;
		// assert(current_run<R);
		//number of c runs before the current run
		uint_t rk = bwt_heads.rank(current_run,c);
		//number of c before i in the current run
		uint_t tail = (bwt_heads[current_run]==c)*dist;
		// in this case, either there are no c before position i
		// or the current run is the first of a certin character
		if(rk==0) return tail;

		return letter_bv[c].select1(rk-1)+1+tail;
	}

	/*
	uint_t rank_(uint_t i, char c, uint_t B){
		// get current run
		uint_t last_block = main_bv.rank1(i);
		uint_t current_run = last_block*B;
		// get first position of the previous block
		uint_t pos = 0;
		if( last_block>0 ){ pos = main_bv.select1(last_block-1)+1; }
		// get distance between i and previous block
		// assert(pos <= i);
		uint_t dist = i-pos;
		//otherwise, scan at most B runs
		while(pos < i){
			// get current run length
			//pos += run_length(current_run);
			pos += run_at(current_run);
			current_run++;
			// update the distance until we get to the
			// correct run
			if(pos<=i) dist = i-pos;
		}
		// get the correct run counter
		if(pos>i) current_run--;
		// assert(current_run<R);
		//number of c runs before the current run
		uint_t rk = bwt_heads.rank(current_run,c);
		//number of c before i in the current run
		uint_t tail = (bwt_heads[current_run]==c)*dist;
		// in this case, either there are no c before position i
		// or the current run is the first of a certin character
		if(rk==0) return tail;

		return letter_bv[c].select1(rk-1)+1+tail;
	}*/

	/*
	 * position of i-th character c. i starts from 0!
	 */
	uint_t select(uint_t i, char c, uint_t B){
		// number of 1s before i
		uint_t j = letter_bv[c].rank1(i);
		//starting position of i-th c inside its run
		//assert(j==0 || i >= runs_per_letter[c].select(j-1) + 1);
		uint_t before = (j==0 ? i : i - (letter_bv[c].select1(j-1) + 1));
		//position in run_heads
		uint_t r = bwt_heads.select(j+1,c);
		//k = number of bits before position of interest in the main string
		//here, k is initialized looking at the sampled runs
		//assert(r/B==0 || r/B-1<runs.number_of_1());
		uint_t k = (r/B==0?0 : main_bv.select1(r/B-1)+1);
		//now add remaining run lengths to k
		for( uint_t t = (r/B)*B; t<r; ++t ){

			k += run_at(t);

		}
		return k + before;
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	uint_t serialize(std::ostream& out){

		uint_t w_bytes = 0;

		out.write((char*)&BWTlength,sizeof(BWTlength));
		out.write((char*)&R,sizeof(R));
		out.write((char*)&B,sizeof(B));
		out.write((char*)C.data(),128*sizeof(uint_t));

		w_bytes += sizeof(BWTlength) + sizeof(R) + sizeof(B) + 128*sizeof(uint_t);

		if(BWTlength==0) return w_bytes;

		w_bytes += main_bv.serialize(out);

		std::vector<int> selChar; int nChar = 0;
		for(int i=0; i<128; ++i){ if( letter_bv[i].size() > 0 ){ selChar.push_back(i); } }
		nChar = selChar.size();
		//std::cout << "Select char: " << nChar << "\n";
		out.write((char*)&nChar,sizeof(nChar));
		out.write((char*)selChar.data(),selChar.size()*sizeof(int));

		w_bytes += sizeof(nChar) + selChar.size()*sizeof(int);

		for(uint_t i=0;i<128;++i){
			if( letter_bv[i].size() > 0 ) w_bytes += letter_bv[i].serialize(out);
		}

		w_bytes += bwt_heads.serialize(out);

		return w_bytes;
	}

	/*
	uint_t letter_size(char c)
	{
		return letter_bv[c].size();
	}*/

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&BWTlength,sizeof(BWTlength));
		in.read((char*)&R,sizeof(R));
		in.read((char*)&B,sizeof(B));
		// load C vector
		C = std::vector<uint_t>(128);
		// C_p = std::vector<bool>(128,0);
		in.read((char*)C.data(),128*sizeof(uint_t));
		// load main bitvector
		main_bv.load(in);
		// load letter bitvectors
		int nChar = 0;
		in.read((char*)&nChar,sizeof(nChar));
		std::vector<int> selChar; selChar.resize(nChar);
		in.read((char*)selChar.data(),selChar.size()*sizeof(int));
		letter_bv = std::vector<sd_vector>(128);
		for(int j=0; j<selChar.size(); ++j)
			{ letter_bv[selChar[j]].load(in); /*C_p[selChar[j]] = 1;*/ }
		// load BWT heads
		bwt_heads.load(in);

	}

private:
	// heads and lengths vectors
	std::vector<char> heads;
	std::vector<uint_t> lens;
	uint_t BWTlength = 0;
	// main bitvector for all characters with support
	// for rank and select queries
	sd_vector main_bv;
	// vector containing one bitvector for each char
	// with support for rank and select queries
	std::vector< sd_vector > letter_bv;
	// number of runs
	uint_t R;
	// block size
	uint_t B;
};



#endif
