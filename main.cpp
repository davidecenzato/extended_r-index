#include <string>
#include <iostream>
#include <chrono>
#include <getopt.h>

// algorithms for computing different BWT variants
#include "r_index.hpp"
// fasta reader function
#include "IOfunc.hpp"

// struct containing command line parameters and other globals
struct args {
  uint_t B = 2;
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
  bool pocc = false;
};

// function that prints the instructions for using the tool
void print_help(char** argv) {
  std::cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << std::endl;
  std::cout << "  Options: " << std::endl
        << "\t-c \tconstruct and store ebwt r-index, def. False" << std::endl
        << "\t-q \tcompute count/locate queries ( 0 (count) | 1 (cout print no. occ.) | 2 (locate) | 3 (locate print occ.) ), def. -1" << std::endl
        << "\t-b B\tbitvector block size, def. 2" << std::endl
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
  if(!arg.build && (arg.query < 0 || arg.query > 3 ) ){ std::cerr << "Error! select a correct mode (either -c | -q 0 | -q 1 | -q 2 | -q 3).\n";  }
}

int main(int argc, char** argv)
{
  // translate command line arguments
  args arg;
  parseArgs(argc, argv, arg);
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

    size_t query_time = 0;

    auto t3 = std::chrono::high_resolution_clock::now();

    if(arg.query < 2){
      FILE * nocc, * ptime;
      arg.pocc = (arg.query == 1);
      std::cout << "Computing count queries..." << std::endl;
      if(arg.query == 1)
      {
        std::string output_file   = arg.patname + ".noccEBWT";
        std::string output_file2  = arg.patname + ".timeEBWT";
        // open output file
        nocc  = fopen(output_file.c_str(), "w+");
        ptime = fopen(output_file2.c_str(),"w+");
      
  	    //extract patterns from file and search them in the index
        //if(arg.pocc){
    		for(int64_t i=0; i<noSeq; ++i){

    			perc = (100*i)/noSeq;
    			if( perc > last_perc ){
    				std::cout << perc << "% done ..." << std::endl;
    				last_perc = perc;
    		  }

    			getline(ifs, pattern);
    			getline(ifs, pattern);

          auto before = std::chrono::high_resolution_clock::now();
    			auto rn = idx.count(pattern);
          uint_t curr_occ = rn.second>=rn.first ? (rn.second-rn.first)+1 : 0;
          auto after = std::chrono::high_resolution_clock::now();

          fwrite(&curr_occ,4,1,nocc);
          std::chrono::duration<double, std::milli> patt_time = after - before;
          float dur_patt = patt_time.count();
          fwrite(&dur_patt,sizeof(float),1,ptime);
          occ_tot += curr_occ;
    		}
        // close output files
        fclose(nocc); fclose(ptime);
      }
      else
      {
        for(int64_t i=0; i<noSeq; ++i){

          perc = (100*i)/noSeq;
          if( perc > last_perc ){
            std::cout << perc << "% done ..." << std::endl;
            last_perc = perc;
          }

          getline(ifs, pattern);
          getline(ifs, pattern);

          auto before = std::chrono::high_resolution_clock::now();
          auto rn = idx.count(pattern);
          auto after = std::chrono::high_resolution_clock::now();
          occ_tot += rn.second>=rn.first ? (rn.second-rn.first)+1 : 0;
          query_time += std::chrono::duration_cast<std::chrono::nanoseconds>(after - before).count();
        }
      }

  		double occ_avg = (double)occ_tot / noSeq;

      STAT[0] = occ_tot; STAT[1] = occ_avg;

  		std::cout << std::endl << occ_avg << " average occurrences per pattern" << std::endl;
    }
    else if(arg.query > 1)
    {
      std::cout << "Computing locate queries..." << std::endl;
      // initialize output file
      FILE * occ;
      if(arg.query==3)
      {
        std::string output_file  = arg.filename + ".occ";
        // open output file
        occ = fopen(output_file.c_str(),"w+");

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

          if(OCC.size() > 0) fwrite(&OCC[0],5,OCC.size(),occ);
      		
      		occ_tot += OCC.size();
    		}
        // close output file
        if(arg.query==3) fclose(occ);
      }
      else
      {
        //extract patterns from file and search them in the index
        for(int64_t i=0; i<noSeq; ++i){

          perc = (100*i)/noSeq;
          if( perc > last_perc ){
            std::cout << perc << "% done ..." << std::endl;
            last_perc = perc;
          }

          getline(ifs, pattern);
          getline(ifs, pattern);

          auto before = std::chrono::high_resolution_clock::now();
          auto OCC = idx.locate_all(pattern,arg.first);
          auto after = std::chrono::high_resolution_clock::now();
          query_time += std::chrono::duration_cast<std::chrono::nanoseconds>(after - before).count();
          
          occ_tot += OCC.size();

        }
      }

  		double occ_avg = (double)occ_tot / noSeq;

      STAT[0] = occ_tot; STAT[1] = occ_avg;

  		std::cout << std::endl << occ_avg << " average occurrences per pattern" << std::endl;

  	}

    ifs.close();

    auto t4 = std::chrono::high_resolution_clock::now();

    uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Load time : " << load << " milliseconds" << std::endl;

    uint64_t search = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
    std::cout << "number of patterns n = " << noSeq << std::endl;
    std::cout << "total number of occurrences  occ_t = " << occ_tot << std::endl;

    std::cout << "Total time : " << (double)query_time/1000000 << " milliseconds" << std::endl;
    std::cout << "Search time : " << ((double)query_time/1000000)/noSeq << " milliseconds/pattern (total: " << noSeq << " patterns)" << std::endl;
    std::cout << "Search time : " << ((double)query_time/1000000)/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << std::endl;
    
    // store some statistics
    STAT[2] = (double)query_time/1000000; STAT[3] = ((double)query_time/1000000)/noSeq; STAT[4] = ((double)query_time/1000000)/occ_tot;
    std::string stat_file  = arg.filename + ".stats";
    FILE * stat = fopen(stat_file.c_str(),"w");
    fwrite(&STAT[0],sizeof(double),5,stat);
    fclose(stat);
  }

  return 0;
}