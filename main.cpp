#include <string>
#include <iostream>
#include <chrono>
#include <getopt.h>

// algorithms for computing different BWT variants
#include "r_index.hpp"
// memory counter
#include "external/malloc_count/malloc_count.h"
// fasta reader function
#include "IOfunc.hpp"

// struct containing command line parameters and other globals
struct args {
  uint_t B = 1;
  bool pfpebwt = false;
  std::string filename = "";
  std::string outname = "";
  std::string patname = "";
  bool read_from_stream = false;
  bool verbose = false; // verbosity level
  bool build = false;
  int query = -1;
  bool check = false; // debug only
  bool first = false;
};

// function that prints the instructions for using the tool
void print_help(char** argv) {
  std::cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << std::endl;
  std::cout << "  Options: " << std::endl
        << "\t-c \tconstruct and store ebwt r-index, def. False" << std::endl
        << "\t-q \tcompute count/locate queries ( 0 (count) | 1 (locate) | 2 (store occurrences) ), def. -1" << std::endl
        << "\t-b B\tbitvector block size, def. 1" << std::endl
        //<< "\t-i \tread pfpebwt files, def. false " << std::endl
        //<< "\t-s \tread input data from stream, def. False " << std::endl
        << "\t-f \tsampled first rotations, def. False " << std::endl
        << "\t-v \tset verbose mode, def. False " << std::endl
        << "\t-p P\tpattern file path, def. <input filename.pat> " << std::endl
        << "\t-o O\tbasename for the output files, def. <input filename>" << std::endl
        << "\t-d \tcheck locate output (debug only)" << std::endl;

  exit(-1);
}
// function for parsing the input arguments
void parseArgs( int argc, char** argv, args& arg ) {
  int c;
  extern int optind;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");
 
  std::string sarg;
  while ((c = getopt( argc, argv, "b:o:q:p:vcsihdf") ) != -1) {
    switch(c) {
      case 'c':
        arg.build = true; break;
        // build index
      case 'v':
        arg.verbose = true; break;
        // verbose mode
      case 's':
        arg.read_from_stream = true; break;
        // stream mode
      case 'i':
        arg.pfpebwt = true; break;
        // bigbwt int width
      case 'q':
        //sarg.assign( optarg );
        arg.query = atoi( optarg ); break;
        // store query type
      case 'p':
        sarg.assign( optarg );
        arg.patname.assign( sarg ); break;
        // store the pattern file path
      case 'b':
        //sarg.assign( optarg );
        arg.B = atoi( optarg ); break;
        // store the bitvector block size
      case 'o':
        sarg.assign( optarg );
        arg.outname.assign( sarg ); break;
        // store the output files path
      case 'd':
        arg.check = true; break;
        // check locate output
      case 'f':
        arg.first = true; break;
        // check locate output
      case 'h':
        print_help(argv); exit(-1);
        // fall through
      default:
        std::cout << "Unknown option. Use -h for help." << std::endl;
        exit(-1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind+1) {
    arg.filename.assign( argv[optind] );
  }
  else {
    std::cout << "Invalid number of arguments" << std::endl;
    print_help(argv);
  }
  // set output files basename
  if(arg.outname == "") arg.outname = arg.filename;
  // set pattern file path
  if(arg.patname == "") arg.patname = arg.filename+".pat";
  // check mode
  if(!arg.build && (arg.query < 0 || arg.query > 2 ) ){ std::cerr << "Error! select a correct mode (either -c | -q 0 | -q 1 | -q 2).\n";  }
}

int main(int argc, char** argv)
{
  // translate command line arguments
  args arg;
  parseArgs(argc, argv, arg);
  int c = 0;
  assert(c != 0);
  // compute and store the r-Index of the eBWT
  if(arg.build){
    if( arg.verbose ){
      std::cout << "Computing the eBWT r-index of: " << arg.filename 
      << " Main bitvector blocksize selected: " << arg.B << "\n";
      /*if( arg.pfpebwt ){ std::cout << "Reading pfpebwt files\n"; }
      if( arg.pfpebwt || arg.read_from_stream )
      {
        std::cout << "Reading input files from stream\n";
      }*/
    }
    // compute and store the ebwt r-index
    r_index(arg.filename,arg.B,arg.read_from_stream,1,arg.verbose,arg.first);
  }
  else if(!arg.check){

    // load r-index data structures
    std::string input_file  = arg.filename + ".eri";
    // open stream
    std::ifstream in(input_file);

    // start test
    auto t1 = std::chrono::high_resolution_clock::now();

    r_index idx = r_index();
    // load
    idx.load(in);

    auto t2 = std::chrono::high_resolution_clock::now();

    in.close();

    std::cout << "Searching patterns in file: " << arg.patname << std::endl;
    std::ifstream ifs(arg.patname);

    int64_t noSeq = 0; std::string line;
    while (getline(ifs, line)){ noSeq++; }
    noSeq /= 2;
    ifs.clear();
    ifs.seekg(0, std::ios::beg);

    int64_t perc = 0, last_perc = 0;
    int64_t occ_tot=0;
    std::string pattern = std::string();

    // initialize stats vector
    std::vector<double> STAT(5,0);

    auto t3 = std::chrono::high_resolution_clock::now();

    if(arg.query==0){
      std::cout << "Computing count queries..." << std::endl;
  	  //extract patterns from file and search them in the index
  		for(int64_t i=0; i<noSeq; ++i){

  			perc = (100*i)/noSeq;
  			if( perc > last_perc ){
  				std::cout << perc << "% done ..." << std::endl;
  				last_perc = perc;
  			}

  			getline(ifs, pattern);
  			getline(ifs, pattern);

  			auto rn = idx.count(pattern);
  			occ_tot += rn.second>=rn.first ? (rn.second-rn.first)+1 : 0;
  		}

  		double occ_avg = (double)occ_tot / noSeq;

      STAT[0] = occ_tot; STAT[1] = occ_avg;

  		std::cout << std::endl << occ_avg << " average occurrences per pattern" << std::endl;
    }
    else if(arg.query>=1)
    {
      std::cout << "Computing locate queries..." << std::endl;
      // initialize output file
      FILE * occ;
      if(arg.query==2)
      {
      std::string output_file  = arg.filename + ".occ";
      // open output file
      occ = fopen(output_file.c_str(),"w+");
      }

  		//extract patterns from file and search them in the index
  		for(int64_t i=0; i<noSeq; ++i){

    		perc = (100*i)/noSeq;
    		if( perc > last_perc ){
    			std::cout << perc << "% done ..." << std::endl;
    			last_perc = perc;
    		}

    		getline(ifs, pattern);
    		getline(ifs, pattern);

    		auto OCC = idx.locate_all(pattern,arg.first);

        if(arg.query==2 && OCC.size() > 0) fwrite(&OCC[0],5,OCC.size(),occ);
    		
    		occ_tot += OCC.size();

  		}

  		double occ_avg = (double)occ_tot / noSeq;

      STAT[0] = occ_tot; STAT[1] = occ_avg;

  		std::cout << std::endl << occ_avg << " average occurrences per pattern" << std::endl;

  		if(arg.query==2) fclose(occ);
  	}

    ifs.close();

    auto t4 = std::chrono::high_resolution_clock::now();

    uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Load time : " << load << " milliseconds" << std::endl;

    uint64_t search = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
    std::cout << "number of patterns n = " << noSeq << std::endl;
    //cout << "pattern length m = " << m << endl;
    std::cout << "total number of occurrences  occ_t = " << occ_tot << std::endl;

    std::cout << "Total time : " << search << " milliseconds" << std::endl;
    std::cout << "Search time : " << (double)search/noSeq << " milliseconds/pattern (total: " << noSeq << " patterns)" << std::endl;
    std::cout << "Search time : " << (double)search/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << std::endl;

    STAT[2] = search; STAT[3] = (double)search/noSeq ; STAT[4] = (double)search/occ_tot;

    std::string stat_file  = arg.filename + ".stats";
    FILE * stat = fopen(stat_file.c_str(),"w");
    fwrite(&STAT[0],sizeof(double),5,stat);
    fclose(stat);

  }
  else
  {
  	std::cout << "testing locate output correctness...\n";
  	 // load r-index data structures
    std::string input_file  = arg.filename + ".eri";
    std::cout << input_file << "\n";
    // open stream
    std::ifstream in(input_file);
    // init empty index
    r_index idx = r_index();
	  // load r-index
    idx.load(in);

    // load input text
    std::vector<unsigned char> Text;
    uint_t n=0, ns=0;
    std::vector<uint_t> onset = {0};
    // initialize gap encoded bit-vector 
    sdsl::sd_vector<> b;
    sdsl::sd_vector<>::rank_1_type r_s;
    sdsl::sd_vector<>::select_1_type s_s;
    // load fasta
    load_fasta(arg.filename.c_str(),Text,onset,n,ns,false,idx.getBWTlen());
    // construct 
    sdsl::sd_vector_builder builder(n+1,onset.size());
    for(auto idx: onset){ builder.set(idx); }
    b = sdsl::sd_vector<>(builder);
    // clear onset vector
    onset.clear();
    r_s = sdsl::sd_vector<>::rank_1_type(&b); 
    s_s = sdsl::sd_vector<>::select_1_type(&b);
    std::cout << "Locating patterns ... " << std::endl;
	  std::ifstream ifs(arg.patname);

    int64_t noSeq = 0; std::string line;
    int64_t perc = 0, last_perc = 0;
    while (getline(ifs, line)){ noSeq++; }
    noSeq /= 2;
    ifs.clear();
    ifs.seekg(0, std::ios::beg);

    std::string pattern = std::string();

    std::cout << "Number of patterns: " << noSeq << std::endl;

  	//extract patterns from file and search them in the index
  	for(int64_t i=0; i<noSeq; ++i){

  		perc = (100*i)/noSeq;
  		if( perc > last_perc ){
  			std::cout << perc << "% done ..." << std::endl;
  			last_perc = perc;
  		}

  		getline(ifs, pattern);
  		getline(ifs, pattern);

  		auto OCC = idx.locate_all(pattern,arg.first);
      uint_t last_occ;
      if(arg.first){ last_occ = idx.Phi_first(OCC[OCC.size()-1]); }
      else{ last_occ = idx.Phi(OCC[OCC.size()-1]); }

      // search for duplicated
      std::sort(OCC.begin(), OCC.end());
      auto it = std::unique(OCC.begin(), OCC.end());
      bool wasUnique = (it == OCC.end());
      if(!wasUnique){ std::cerr << "Duplicate indexes detected! pattern no: " << i << std::endl; exit(1); }
      // last_occ must not contain the pattern
      OCC.push_back(last_occ);

      // check located indexes correctness
  		for(int64_t i=0; i<OCC.size(); ++i){
        bool flag = false; 
  			int64_t next = s_s(r_s(OCC[i]+1)+1);
  			if( ( next - OCC[i] ) < pattern.size()  )
  			{
  				int64_t j = 0;
  				for(; j<pattern.size(); ++j)
  				{
  					if( next == OCC[i]+j ){ break; }
  					if( pattern[j] != Text[OCC[i]+j] )
            {
              if(i < OCC.size()-1){ std::cerr << OCC[i] << " " << pattern << " Err_1!\n"; exit(1); }
              else{ flag = true; break; }
            }
  				}
  				//j--;
          if( !flag ){
    				int64_t curr = s_s(r_s(OCC[i]+1));
    				for(; j<pattern.size(); ++j)
    				{
    					if(pattern[j] != Text[curr++] )
              {
                if(i < OCC.size()-1){ std::cerr << OCC[i] << " " << pattern << " Err_2!\n"; exit(1); }
                else{ flag = true; break; }
              }
    				}
          }
  			}
  			else
  			{
  				for(int64_t j=0; j<pattern.size(); ++j)
  				{
  					if( pattern[j] != Text[OCC[i]+j] )
              {
                if(i < OCC.size()-1){ std::cerr << OCC[i] << " " << pattern << " Err_3!\n"; exit(1); }
                else{ flag = true; break; }
              }
  				}
  			} 
        // check last occurrence
        if( i==OCC.size()-1 && !flag ){ std::cerr << "The interval can be extended...\n"; exit(1); }
  		}
  	}
    std::cout << "100% done ...\nEverything's fine!" << std::endl;
  }

  return 0;
}