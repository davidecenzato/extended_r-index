#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path

Description = """
Tool to build the extended r-index of string collections.
"""

dirname = os.path.dirname(os.path.abspath(__file__))
extrindex_exe     =  os.path.join(dirname, "build/er-index")
extrindex64_exe   =  os.path.join(dirname, "build/er-index64")
parseNT_exe       =  os.path.join(dirname, "build/circpfpNT.x")
parsebwtNT_exe    =  os.path.join(dirname, "build/parsebwtNT.x")
bebwtNT_exe       =  os.path.join(dirname, "build/bebwtNT.x")
bebwtNTp64_exe    =  os.path.join(dirname, "build/bebwtNTp64.x")
bebwtNTd64_exe    =  os.path.join(dirname, "build/bebwtNTd64.x")
bebwtNT64_exe     =  os.path.join(dirname, "build/bebwtNT64.x")

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input', help='input fasta file name', type=str)
    parser.add_argument('--construct', help='constructs the extended r-index (def. False)', action='store_true')
    parser.add_argument('-w', '--wsize', help='sliding window size for PFP (def. 10)', default=10, type=int)
    parser.add_argument('-p', '--mod', help='hash modulus for PFP (def. 100)', default=100, type=int)
    parser.add_argument('-b', '--B', help='bitvector block size for predecessor queries (def. 2)', default=2, type=int)
    parser.add_argument('--first', help='sample first rotation of each sequence (def. False)', action='store_true')
    #parser.add_argument('-a', '--algo', help='eBWT construction algorithm (def. bigbwt)', default="bigbwt", type=str)
    #parser.add_argument('-t', help='number of helper threads (def. None)', default=0, type=int)
    #parser.add_argument('-n', help='number of different primes (def. 1)', default=1, type=int)
    #parser.add_argument('--query', help='compute count and locate queries (def. False)', action='store_true')
    parser.add_argument('--pfile', help='pattern file path (def. <input filename.pat>)', default="", type=str)
    parser.add_argument('--count', help='compute count queries (def. False)', action='store_true')
    parser.add_argument('--locate', help='compute locate queries (def. False)', action='store_true')
    parser.add_argument('--verbose',  help='verbose (def. False)',action='store_true')
    #parser.add_argument('-d',  help='use remainders instead of primes (def. False)',action='store_true')
    #parser.add_argument('--reads', help='process input ad a reads multiset (def. False)', action='store_true')
    #parser.add_argument('--period', help='remove periodic sequences (def. False)', action='store_true')
    args = parser.parse_args()

    logfile_name = args.input + ".log"
    # get main bigbwt directory
    args.extrindex_dir = os.path.split(sys.argv[0])[0]
    print("Sending logging messages to file:", logfile_name)
    with open(logfile_name,"a") as logfile:

        #if(args.algo == "bigbwt"):
        if( args.construct ):
            # ---------- Parsing of the input file
            start0 = start = time.time()
            '''
            if args.reads:
                # Input is a short sequences multiset
                if args.d:
                    # Use different remainders
                    if args.t>0:
                        command = "{exe} {file} -w {wsize} -p {modulus} -t {th} -n {wnumb}".format(
                                exe = os.path.join(args.extrindex_dir,parseReadsD_exe),
                                wsize=args.wsize, modulus = args.mod, th=args.t, wnumb=args.n, file=args.input)
                    else:
                        command = "{exe} {file} -w {wsize} -p {modulus} -n {wnumb}".format(
                                exe = os.path.join(args.extrindex_dir,parseReadsDNT_exe),
                                wsize=args.wsize, modulus = args.mod, wnumb = args.n, file=args.input)
                else:
                    # Use different primes
                    if args.t>0:
                        command = "{exe} {file} -w {wsize} -p {modulus} -t {th} -n {wnumb}".format(
                                exe = os.path.join(args.extrindex_dir,parseReads_exe),
                                wsize=args.wsize, modulus = args.mod, th=args.t, wnumb=args.n, file=args.input)
                    else:
                        command = "{exe} {file} -w {wsize} -p {modulus} -n {wnumb}".format(
                                exe = os.path.join(args.extrindex_dir,parseReadsNT_exe),
                                wsize=args.wsize, modulus = args.mod, wnumb = args.n, file=args.input)
            else:
            '''
                # Input is a long sequences multiset
            '''
            if args.t>0:
                command = "{exe} {file} -w {wsize} -p {modulus} -t {th}".format(
                        exe = os.path.join(args.extrindex_dir,parse_exe),
                        wsize=args.wsize, modulus = args.mod, th=args.t, file=args.input)
            else:
            '''
            command = "{exe} {file} -w {wsize} -p {modulus}".format(
                    exe = os.path.join(args.extrindex_dir,parseNT_exe),
                    wsize=args.wsize, modulus = args.mod, file=args.input)

            print("==== Parsing. Command:", command)
            if(execute_command(command,logfile,logfile_name)!=True):
                return
            print("Elapsed time: {0:.4f}".format(time.time()-start))

            # ----------- Compute the inverted list of parse's ebwt
            start = time.time()
            print("==== Computing Inverted list of parse's eBWT.")
            parse_size = os.path.getsize(args.input+".eparse")/4
            print("Parse contains " + str(parse_size) + " words.")
            if(parse_size >= (2**32-1)):
                print("IL creation running in 64 bit mode")
                command = "{exe} {file} -w {wsize}".format(
                        exe = os.path.join(args.extrindex_dir,parsebwtNT64_exe), wsize=args.wsize, file=args.input)
            else:
                print("IL creation running in 32 bit mode")
                command = "{exe} {file} -w {wsize}".format(
                         exe = os.path.join(args.extrindex_dir,parsebwtNT_exe), wsize=args.wsize, file=args.input)

            print("Command:", command)
            if(execute_command(command,logfile,logfile_name)!=True):
                return
            print("Elapsed time: {0:.4f}".format(time.time()-start))

            # ----------- Computing the eBWT of the text
            start = time.time()
            print("==== Computing the eBWT and samples of the text.")
            dict_size = os.path.getsize(args.input+".edict")
            print("Dictionary contains " + str(dict_size) + " characters.")
            if(dict_size >=  (2**31-1)):
                print("Dict SA running in 64 bit mode")
                if(parse_size >= (2**32-1)):
                    command = "{exe} {file} -w {wsize}".format(
                            exe = os.path.join(args.extrindex_dir,bebwtNT64_exe), wsize=args.wsize, file=args.input)
                else:
                    command = "{exe} {file} -w {wsize}".format(
                            exe = os.path.join(args.extrindex_dir,bebwtNTd64_exe), wsize=args.wsize, file=args.input)
            else:
                print("Dict SA running in 32 bit mode")
                if(parse_size >= (2**32-1)):
                    command = "{exe} {file} -w {wsize}".format(
                            exe = os.path.join(args.extrindex_dir,bebwtNTp64_exe), wsize=args.wsize, file=args.input)
                else:
                    command = "{exe} {file} -w {wsize}".format(
                            exe = os.path.join(args.extrindex_dir,bebwtNT_exe), wsize=args.wsize, file=args.input)
            # output the eBWT in rle format
            command += " -r"
            # output the GCA-samples
            command += " -s"
            # sample the first rotation of each sequence
            if(args.first): command += " -f"
            print("==== Computing the eBWT and the GCA-samples of the input. Command:", command)
            if(execute_command(command,logfile,logfile_name)!=True):
                return
            # print total time
            print("Elapsed time: {0:.4f}".format(time.time()-start));
            #print("Total construction time: {0:.4f}".format(time.time()-start0))
            # delete auxiliary files
            print("Deleting auxiliary files")
            command = "rm -f {file}.eparse_old {file}.offset_old {file}.eparse {file}.edict {file}.offset {file}.eocc {file}.fchar {file}.start {file}.sdsl".format(file=args.input)
            if(execute_command(command,logfile,logfile_name)!=True):
                return

            start = time.time()
            ## construct extended r-index
            input_size = os.path.getsize(args.input)
            f = open(args.input+".mode", "w")
            ## construct command
            if(input_size >=  (2**32-1)):
                command = "{exe} {file} -c -b {bsize}".format(
                exe = os.path.join(args.extrindex_dir,extrindex64_exe),
                bsize=args.B, file=args.input)
                f.write("64")
            else:
                command = "{exe} {file} -c -b {bsize}".format(
                exe = os.path.join(args.extrindex_dir,extrindex_exe),
                bsize=args.B, file=args.input)
                f.write("32")
            f.close()
            # sample the first rotation of each sequence
            if(args.first): command += " -f"
            # execute command
            print("==== Computing the extended r-index of the input. Command:", command)
            if(execute_command(command,logfile,logfile_name)!=True):
                return
            print("Elapsed time: {0:.4f}".format(time.time()-start));
            print("Total construction time: {0:.4f}".format(time.time()-start0))

            index_size = os.path.getsize(args.input+".eri")
            print("Original input size: " + str(input_size) + " bytes" )
            print("Extended r-index size: " + str(index_size) + " bytes" )

        ## queries
        print("qua")
        if( args.count ):
            f = open(args.input+".mode", "r")
            mode = int(f.read())
            f.close()
            if(mode > 32):
                command = "{exe} {file} -q 0 -p {pfile} ".format(
                exe = os.path.join(args.extrindex_dir,extrindex64_exe),
                file=args.input, pfile=args.pfile)
            else:
                command = "{exe} {file} -q 0 -p {pfile} ".format(
                exe = os.path.join(args.extrindex_dir,extrindex_exe),
                file=args.input, pfile=args.pfile)
            if(args.first): command += " -f"
            print("==== Computing count queries. Command:", command)
            subprocess.check_call( command.split() )

        if( args.locate ):
            f = open(args.input+".mode", "r")
            mode = int(f.read())
            f.close()
            if(mode > 32):
                command = "{exe} {file} -q 2 -p {pfile} ".format(
                exe = os.path.join(args.extrindex_dir,extrindex64_exe),
                file=args.input, pfile=args.pfile)
            else:
                command = "{exe} {file} -q 2 -p {pfile} ".format(
                exe = os.path.join(args.extrindex_dir,extrindex_exe),
                file=args.input, pfile=args.pfile)
            if(args.first): command += " -f"
            print("==== Computing locate queries. Command:", command)
            subprocess.check_call( command.split() )


# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name,env=None):
  try:
    #subprocess.run(command.split(),stdout=logfile,stderr=logfile,check=True,env=env)
    subprocess.check_call(command.split(),stdout=logfile,stderr=logfile,env=env)
  except subprocess.CalledProcessError:
    print("Error executing command line:")
    print("\t"+ command)
    print("Check log file: " + logfile_name)
    return False
  return True
  

if __name__ == '__main__':
    main()
