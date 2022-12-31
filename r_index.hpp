#ifndef R_INDEX_S_H_
#define R_INDEX_S_H_

#include <string>
#include <vector>
#include <iostream>
#include <cassert>

#include <sdsl/wavelet_trees.hpp>
#include "rle_ebwt.hpp"
#include "pred_ebwt.hpp"

typedef rle_ebwt rle_t;
typedef pred_ebwt pred_t;
typedef std::pair<uint_t,uint_t> range_t;

// define r index class
class r_index{

public:
	// empty constructor
	r_index(){}
	// constructor
	r_index(std::string input, uint_t bsize = 1, bool stream = 0, bool pfpebwt = 0, bool verbose = 0, bool first = 0){
		// get int size
		int isize = sizeof(uint_t);
		if( pfpebwt ){ isize = 5; }

		std::cout << "(1/3) Compute the RLE eBWT data structure\n";
		// input files
		std::string heads = input + ".head";
		std::string lens = input + ".len";
		// set block size
		B = bsize;
		if(!stream && !pfpebwt){
			// run length encoded eBWT
			bwt = rle_ebwt(heads, lens, B,verbose);
		}
		else{
			// open streams
			std::ifstream head_s(heads);
			std::ifstream len_s(lens);
			// run length encoded eBWT
			bwt = rle_ebwt(head_s, len_s, heads, B, isize,verbose);
		}

		std::cout << "(2/3) Compute the predecessor search data structure\n";
		// input files
		std::string s_samples = input + ".ssam";
		std::string e_samples = input + ".esam";
		std::string st_pos = input + ".spos";
		
		if(!stream && !pfpebwt){
			// construct predecessor data structure for the eBWT
			phi = pred_ebwt(s_samples, e_samples, st_pos, verbose);
			//phi.construct_rank_select_dt();
		}
		else{
			// open streams
			std::ifstream s_samples_s(s_samples);
			std::ifstream e_samples_s(e_samples);
			std::ifstream st_pos_s(st_pos);
			// construct predecessor data structure for the eBWT
			phi = pred_ebwt(s_samples_s, e_samples_s, st_pos_s, bwt.size(), isize, verbose, first);
			//phi.construct_rank_select_dt();
		}

        std::cout << "(3/3) Serialize the eBWT r-index data structure\n";
		std::string path = input.append(".eri");
		std::ofstream out(path);

		uint_t space = serialize(out);
		if(verbose) std::cout << "TOT space: " << space << " Bytes" << std::endl << std::endl;

		out.close();

	}

   /*
	* Return eBWT range of pattern P
	*/
	range_t count(std::string &P){

		range_t range = {0,bwt.size()-1};
		uint_t m = P.size();

		for(int i=0; i<m and range.second >= range.first; ++i){
			// get new range
			range = LF(range,P[m-i-1]);
		}

		return range;
	}

	/*
	 * return the new range computed applying a LF step to range rn using character c
	 */
	range_t LF(range_t rn, char c){

		//if character does not appear in the text, return empty pair
		if((c==127 and bwt.C[c]==bwt.size()) || bwt.C[c]>=bwt.C[c+1]){ return {1,0}; }
		//number of c before the interval
		uint_t c_before = bwt.rank(rn.first,c,B);
		//number of c inside the interval rn
		uint_t c_inside = bwt.rank(rn.second+1,c,B) - c_before;
		//if there are no c in the interval, return empty range
		if(c_inside==0) return {1,0};

		uint_t l = bwt.C[c] + c_before;

		return {l,l+c_inside-1};
	}
	

	/*
	 * return the predecessor of i in Conjugate array order
	 */
	uint_t Phi(uint_t i){

		//jr is the rank of the predecessor of i (circular)
		//auto pred_query = phi.circular_rank_predecessor_triple(i);
		auto pred_query = phi.circular_rank_predecessor_tuple(i);
		uint_t jr = std::get<0>(pred_query);

		// the actual predecessor
		//uint_t j = phi.select(jr);
		uint_t j = std::get<1>(pred_query);

		// compute distance between the two indices
		uint_t delta = 0;
		if( j <= i ){ delta = i-j; }
		else
		{ 
			delta += i-std::get<2>(pred_query);
			delta += std::get<3>(pred_query)-j+1;
		}

		// sample at the end of previous run
		uint_t prev_sample = phi.sample_last( phi.f_to_r(jr)-1 );
		// get starting position of the next sequence
		uint_t next = phi.next_start_pos(prev_sample);

		if( (prev_sample + delta) < next ){ return prev_sample + delta; }
		else{                               return phi.curr_start_pos(prev_sample) + ((prev_sample + delta)%next); }
	}

	/*
	 * return the predecessor of i in Conjugate array order
	 */
	uint_t Phi_first(uint_t i){

		// jr is the rank of the predecessor of i (circular)
		uint_t jr = phi.circular_rank_predecessor_first(i);

		// the actual predecessor
		uint_t j = phi.select(jr);

		// compute distance between the two indices
		// assert(i >= j);
		uint_t delta = i-j;

		// sample at the end of previous run
		uint_t prev_sample = phi.sample_last( phi.f_to_r(jr)-1 );
		// get starting position of the next sequence
		uint_t next = phi.next_start_pos(prev_sample);
		////std::cout << "next: " << next << "\n";

		if( (prev_sample + delta) < next ){ return prev_sample + delta; }
		else{                               return phi.curr_start_pos(prev_sample) + ((prev_sample + delta)%next); }
	}

	/*
	 * return the conjugate array interval of pattern P + the last sample of the interval
	 */
	std::pair<range_t, uint_t> count_and_get_occ(std::string &P){

		uint_t k = 0, ks = 0;
		char c;
		range_t range = {0,bwt.size()-1};
		k = phi.sample_last(bwt.nrun()-1);
		ks = phi.curr_start_pos(k);
		range_t range1;

		uint_t m = P.size();

		for(uint_t i=0;i<m and range.second>=range.first;++i){
			// current character
			c = P[m-i-1];
			// new range computed with the LF step
			range1 = LF(range,c);
			//if suffix can be left-extended with char
			if(range1.first <= range1.second){
				// compute the last sample of the new interval 
				if(bwt[range.second] == c){
					// last c is at the end of range.
					if( k > ks ){	k--;	}
					else
					{
						k = phi.next_start_pos(k) - 1;
					}
				// else find new sample
				}else{
					// find last c in range (there must be one because range1 is not empty)
					// and get its sample (must be sampled because it is at the end of a run)
					// note: by previous check, bwt[range.second] != c, so we can use argument range.second
					uint_t rnk = bwt.rank(range.second,c,B);
					//this is the rank of the last c
					rnk--;
					//jump to the corresponding BWT position
					uint_t j = bwt.select(rnk,c,B);
					//run of position j
					uint_t run_of_j = bwt.run_of_position(j);
					// get sample
					k = phi.sample_last(run_of_j);
					// get new starting pos
					ks = phi.curr_start_pos(k);
					if( k != ks ){ k--; }
					else
					{
						k = phi.next_start_pos(k)-1;
					}
				}
			}
			range = range1;
		}
		return {range, k};
	}

	/*
	 * locate all occurrences of P and return them in an array
	 * (space consuming if result is big).
	 */
	std::vector<uint_t> locate_all(std::string& P, bool first = 0){

		std::vector<uint_t> OCC;

		std::pair<range_t, uint_t> res = count_and_get_occ(P);

		uint_t L = std::get<0>(res).first;
		uint_t R = std::get<0>(res).second;
		uint_t k = std::get<1>(res);

		//std::cout << L << " : " << R << " - " << k << "\n";

		uint_t n_occ = R>=L ? (R-L)+1 : 0;
		OCC.reserve(n_occ);
		if(n_occ>0){
			// push last sample as first occurrence 
			OCC.push_back(k);
			// for each other occurrence apply Phi step
			for(uint_t i=1; i<n_occ; ++i){
				// compute predecessor gCA value
				if(first){ k = Phi_first(k); }
				else{ k = Phi(k); }
				// add occurrence
				OCC.push_back(k);
			}
		}

		return OCC;
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	uint_t serialize(std::ostream& out){

		uint_t w_bytes = 0;

		out.write((char*)&B,sizeof(B));

		w_bytes += sizeof(B);

		w_bytes += bwt.serialize(out);
		w_bytes += phi.serialize(out);

		return w_bytes;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&B,sizeof(B));

		bwt.load(in);
		phi.load(in);

	}

	uint_t getBWTlen(){
		return bwt.size();
	}

private:
	// run-length encoded eBWT
	rle_t bwt;
	// predecessor data structure eBWT
	pred_t phi;
	// block size
	uint_t B;
};


#endif