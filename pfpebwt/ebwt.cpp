/*
 * Code to build the eBWT of a text using its Prefix-free parse.
 * 
 */

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <assert.h>
#include <errno.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>

extern "C" {
#include "utils.h"
#include "xerrors.h" 
}
#include "common.hpp"
#include "dictionary.hpp"
#include "pfp_parse.hpp"
#include "pfp.hpp"
#include "pfp_ssa.hpp" 

using namespace std;

// struct containing command line parameters and other globals
typedef struct {
   string inputFileName = "";
   int w = 10;
   bool rle = 0;
   bool sample_first = 0;
   bool sample = 0;
} Args;


static void parseArgs(int argc, char** argv, Args *arg ) {
  extern int optind, opterr, optopt;
  extern char *optarg;  
  int c;

  puts("==== Command line:"); 
  for(int i=0;i<argc;i++) 
    printf(" %s",argv[i]); 
  puts("\n");

  while ((c = getopt( argc, argv, "w:rsf") ) != -1) { 
    switch(c) { 
      case 'w':
      arg->w = atoi(optarg); break; 
      case 'r':
      arg->rle = 1; break; 
      case 's':
      arg->sample = true; break;
      case 'f':
      arg->sample_first = true; break;
      case '?':
      puts("Unknown option. Use -h for help.");
      exit(1);
    }
  }

  arg->inputFileName = argv[optind];
}

int main(int argc, char** argv) {
    
    // translate command line parameters 
    Args arg;
    parseArgs(argc, argv, &arg);
    
    // start measuring wall time clock
    time_t start_wc = time(NULL);

    cout << "Loading parse's data structures..." << endl;
    pfp_parse pars(arg.inputFileName);
    
    cout << "Loading parse's data structures took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    start_wc = time(NULL);
    
    cout << "Computing BWT of the dictionary..." << endl;
    dictionary dict(arg.inputFileName,arg.w);
    
    cout << "Building the BWT of the dictionary took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    start_wc = time(NULL);

    if(!arg.sample){ 
      // compute only the eBWT
      cout << "Computing eBWT of the text..." << endl;
      pfp pfp(pars,dict,arg.inputFileName,arg.w,arg.rle);
    }
    else{
      // compute the eBWT + gCA samples
      cout << "Computing the eBWT of the text and the GCA-samples..." << endl;
      pfp_ssa pfp_ssa(pars,dict,arg.inputFileName,arg.w,arg.rle,arg.sample_first);
    }
    
    cout << "Building the eBWT of Text took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    
    return 0;
}