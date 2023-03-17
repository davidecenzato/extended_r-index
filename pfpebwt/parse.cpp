/*
 * Code to build the eBWT of the circular Prefix-free parse.
 * 
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <assert.h>
#include <errno.h>
#include "parse.hpp"

using namespace std;

// struct containing command line parameters and other globals
typedef struct {
   string inputFileName = "";
   int w = 10;
} Args;


static void parseArgs(int argc, char** argv, Args *arg ) {
  extern int optind, opterr, optopt;
  extern char *optarg;  
  int c;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("\n");

  while ((c = getopt( argc, argv, "w:r") ) != -1) {
    switch(c) {
      case 'w':
      arg->w = atoi(optarg); break;
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
    
    cout << "Computing eBWT of the parse..." << endl;
    parse pars(arg.inputFileName);
    
    cout << "Building the eBWT of the parse took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";

    return 0;
}