/*
 * Implementation of the predecessor data structure for the GCA.
 * 
 * This code is adapted from https://github.com/nicolaprezza/r-index.git
 *
 */

#ifndef PRED_EBWT_HPP_
#define PRED_EBWT_HPP_

#include <algorithm>
#include <tuple>
#include <sdsl/int_vector.hpp>
#include "sd_vector.hpp"

/*
 *  give bitsize we need to store x
 */
uint8_t bitsize(uint64_t x){

	if(x==0) return 1;
	return 64 - __builtin_clzll(x);

}

/*
 *  class for sorting the indices of an array
 */
class sort_indices
{
   private:
     uint_t* mparr;
   public:
     sort_indices(uint_t* parr) : mparr(parr) {}
     bool operator()(uint_t i, uint_t j) const { return mparr[i]<mparr[j]; }
};

class pred_ebwt{

public:
	// empty constructor
	pred_ebwt(){}
	/*
 	 *  takes in input the files containin the starting and ending sample
 	 *  and the string offsets and construct predecessor data structure
     */
	pred_ebwt(std::string &s_sample_file, std::string &e_sample_file, std::string &s_pos_file, bool verbose = false){
		// input vectors
		std::vector<uint_t> samples_first_vec;
		std::vector<uint_t> indices;
		std::vector<uint_t> samples_last_vec;
		std::vector<uint_t> onset_vec;
		// read samples first offsets
		read_file(s_pos_file.c_str(),onset_vec);
		uint_t BWT_length = onset_vec[onset_vec.size()-1];
		// construct bit vector of string delimiters
		delim = sd_vector(onset_vec,BWT_length+1);
		// read heads file
		read_file(s_sample_file.c_str(),samples_first_vec);
		uint_t r = samples_first_vec.size();
		indices.reserve(r);
		for(uint_t i=0;i<r;++i){ indices.push_back(i); }
		// sort indices
		std::sort(indices.begin(), indices.end(), sort_indices(&samples_first_vec[0]));
		// compute size necessary to store ending samples and first_to_run data structure
		int log_r = bitsize(uint64_t(r));
		int log_n = bitsize(uint64_t(BWT_length));
		// print some statistics
		if(verbose)
		{
			std::cout << "Value n/(sampled r) = " << double(BWT_length)/r << std::endl;
			std::cout << "Number of bits to store first_to_run vector = " << log_r << std::endl;
			std::cout << "Number of bits to store last samples vector = " << log_n << std::endl;
		}
		// create first_to_run vector and sorted samples vector
		first_to_run = sdsl::int_vector<>(r,0,log_r); 
		for(uint_t i=0;i<r;++i){
			first_to_run[i] = indices[i];
			indices[i] = samples_first_vec[indices[i]];
		}
		// free memory
		samples_first_vec.clear();
		// create compressed bit vector of the sorted first samples
		pred = sd_vector(indices,BWT_length);
		// check sample correctness
		{
			delim.construct_select_ds();
			delim.construct_rank_ds();
			pred.construct_rank_ds();
			uint_t prnk = 0; 
			for(uint_t i=1; i<delim.rank1(delim.size()); ++i){
				uint_t rnk = pred.rank1(delim.select1(i));
				if(prnk == rnk){ std::cerr << "Error in .ssam file, sample missing! Check for string duplicates.\n"; exit(1); }
				prnk = rnk;
			}
		}
		//text positions corresponding to last characters in BWT runs, in BWT order
		read_file(e_sample_file.c_str(),samples_last_vec);
		samples_last = sdsl::int_vector<>(r,0,log_n); 
		// construct last samples data structure
		for(uint_t i=0;i<r;++i){ 
			assert(bitsize(uint64_t(samples_last_vec[i])) <= log_n);
			samples_last[i] = samples_last_vec[i];
		}
		// free memory
		samples_last_vec.clear();
	}

	// 2nd constructor
	pred_ebwt(std::ifstream& s_sample_file, std::ifstream& e_sample_file, std::ifstream& s_pos_file, uint_t BWT_length,
		      int isize, bool verbose = false, bool first = false){
		// input vectors
		std::vector<uint_t> samples_first_vec;
		std::vector<uint_t> indices;
		// set BWT length
		// BWT_length = BWT_length_;
		// construct bit vector of string delimiters
		delim = sd_vector(s_pos_file,BWT_length+1,isize);
		// read heads file
		s_sample_file.seekg(0, std::ios::end);
		uint_t r = s_sample_file.tellg()/isize;
    	s_sample_file.seekg(0, std::ios::beg);
    	// read first samples
		samples_first_vec.resize(r);
		uint64_t currSam = 0;
		for(uint_t i=0;i<r;++i){
			s_sample_file.read(reinterpret_cast<char*>(&currSam), isize);
			samples_first_vec[i] = (uint_t)currSam;
		}
		indices.reserve(r);
		for(uint_t i=0;i<r;++i){ indices.push_back(i); }
		// sort indices
		std::sort(indices.begin(), indices.end(), sort_indices(&samples_first_vec[0]));
		// compute size necessary to store ending samples and first_to_run data structure
		int log_r = bitsize(uint64_t(r));
		int log_n = bitsize(uint64_t(BWT_length));
		// print some statistics
		if(verbose)
		{
			std::cout << "Value n/(sampled r) = " << double(BWT_length)/r << std::endl;
			std::cout << "Number of bits to store first_to_run vector = " << log_r << std::endl;
			std::cout << "Number of bits to store last samples vector = " << log_n << std::endl;
		}
		// create first_to_run vector and sorted samples vector
		first_to_run = sdsl::int_vector<>(r,0,log_r); 
		for(uint_t i=0;i<r;++i){
			first_to_run[i] = indices[i];
			indices[i] = samples_first_vec[indices[i]];
		}
		// free memory
		samples_first_vec.clear();
		// create compressed bit vector of the sorted first samples
		pred = sd_vector(indices,BWT_length);
		// check sample correctness
		if(!first)
		{
			delim.construct_select_ds();
			delim.construct_rank_ds();
			pred.construct_rank_ds();
			uint_t prnk = 0;
			for(uint_t i=1; i<delim.rank1(delim.size()); ++i){
				uint_t rnk = pred.rank1(delim.select1(i));
				if(prnk == rnk){
					std::cerr << "Error in .ssam file, sample missing! Please restart computation using --first flag.\n";
					std::cerr << "Sample missing in string number: " << i << "\n";
					exit(1);
				}
				prnk = rnk;
			}
		}
		else
		{
			delim.construct_select_ds();
			delim.construct_rank_ds();
			for(uint_t i=0; i<delim.rank1(delim.size())-1; ++i){
				if( !pred.at(delim.select1(i)) )
				{ std::cerr << "Error in .ssam file, sample missing in string number: " << i+1 << "\n";
				  exit(1);
				}
			}
		}
		//text positions corresponding to last characters in BWT runs, in BWT order
		samples_last = sdsl::int_vector<>(r,0,log_n); 
		// construct last samples data structure
		uint64_t currLastS = 0;
		for(uint_t i=0;i<r;++i){ 
			// compute last samples vector
			e_sample_file.read(reinterpret_cast<char*>(&currLastS), isize);
			samples_last[i] = (uint_t)currLastS;
		}
		// close stream
		e_sample_file.close();
	}

	/*
 	 *  construct rank select data structures for all bitvectors
 	 */
	/*
	void construct_rank_select_dt(){
		// rank select for main bitvector
		pred.construct_rank_ds();
		pred.construct_select_ds();
		// rank select for main bitvector
		delim.construct_rank_ds();
		delim.construct_select_ds();
	}*/

	/*
 	 *  compute the rank of the circular predecessor of i. Returns a tuple containing
 	 *  the rank of the predecessor, the position of the predecessor, the current string
 	 *  starting point and the starting point of the next string
 	 */
	std::tuple<uint_t,uint_t,uint_t,uint_t> circular_rank_predecessor_tuple(uint_t i){
		// compute number of samples before position i
		uint_t rank = pred.rank1(i+1);
		// if there is no predecessor
		if( rank == 0 ){ 
			// return rank of the first sample occurring
			// at the end of the first string
			uint_t last_pos = delim.select1(1);
			rank = pred.rank1(last_pos);
			return std::make_tuple(rank-1,pred.select1(rank-1),0,last_pos-1);
		}
		else{
			// compute actual predecessor position
			uint_t p_pos = pred.select1(rank-1);
			// compute current string index
			uint_t str_id = delim.rank1(i+1);
			uint_t st_pos = delim.select1(str_id-1);
			// check if the current predecessor is in the correct string
			if( p_pos >= st_pos ){ return std::make_tuple(rank-1,p_pos,0,0); }
			else
			{
				// compute the starting point of the next sequence
				uint_t last_pos = delim.select1(str_id);
				// compute the correct rank
				rank = pred.rank1(last_pos);
				//uint_t nrank = pred.rank1(last_pos);
				return std::make_tuple(rank-1,pred.select1(rank-1),st_pos,last_pos-1);
			}
		}
	} 

	/*
 	 *  compute the rank of the circular predecessor of i. If the first position
 	 *  is always sampled than we can just take the first predecessor.
 	 */
	uint_t circular_rank_predecessor_first(uint_t i){
		// compute number of samples before position i
		return pred.rank1(i+1)-1;
	} 

	/*
		return position of ith bit
	*/
	uint_t select(uint_t i){
		return pred.select1(i);
	}

	/*
		return starting point of the next string
	*/
	uint_t next_start_pos(uint_t i){
		return delim.select1(delim.rank1(i+1));
	}

	/*
		return starting point of the current string
	*/
	uint_t curr_start_pos(uint_t i){
		return delim.select1(delim.rank1(i+1)-1);
	}

	/*
		return no. of runs of the ebwt
	*/
	//uint_t no_runs(){
	//	return r;
	//}

	/*
		return eBWT length
	*/
	//uint_t bwt_len(){
	//	return BWT_length;
	//}

	/*
		return the mapping between a beginning sample
		and the run on which it has been sampled
	*/
	uint_t f_to_r(uint_t i){
		return first_to_run[i];
	}

	/*
		return last sample
	*/
	uint_t sample_last(uint_t i){
		return samples_last[i];
	}

	/*  serialize the structure to the ostream
	 *  \param out	 the ostream
	 */
	uint_t serialize(std::ostream& out){

		uint_t w_bytes = 0;

		//out.write((char*)&BWT_length,sizeof(BWT_length));
		//out.write((char*)&r,         sizeof(r));

		//w_bytes += sizeof(BWT_length) + sizeof(r);


		//if(BWT_length==0) return w_bytes;

		w_bytes += pred.serialize(out);
		w_bytes += delim.serialize(out);
		w_bytes += samples_last.serialize(out);
		w_bytes += first_to_run.serialize(out);

		return w_bytes;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		//in.read((char*)&BWT_length,sizeof(BWT_length));
		//in.read((char*)&r,         sizeof(r));

		pred.load(in);
		delim.load(in);
		samples_last.load(in);
		first_to_run.load(in);
	}

private:
	// the predecessor structure on positions corresponding to first chars in BWT runs
	sd_vector pred, delim;
	// text positions corresponding to last characters in BWT runs, in BWT order
	sdsl::int_vector<> samples_last; 
	// stores the BWT run (0...R-1) corresponding to each position in pred, in text order
	sdsl::int_vector<> first_to_run; 
	// BWT length
	// uint_t BWT_length;
	// no runs
	// uint_t r;

};

#endif