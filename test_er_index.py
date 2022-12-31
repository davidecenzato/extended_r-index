#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, datetime, multiprocessing, signal, re


Description = """
Python script for testing the eBWT r-index
Adapted from a code by Massimiliano Rossi
"""

dirname         = os.path.dirname(os.path.abspath(__file__))

#------------------------Executables--------------------------------------
er_index_exe = "/home/davide/ebwt_r-index/er-index" 
er_index_exe_64 = "/home/davide/ebwt_r-index/er-index64" 
genpattern_exe = "/home/davide/ebwt_r-index/genpattern" 
pfpebwt_exe = "/home/davide/rIndexEBWT/PFP-eBWT/pfpebwt"
r_index_nic_count_exe = "/home/davide/r-indexNicola/r-index/build/ri-count"
r_index_nic_locate_exe = "/home/davide/r-indexNicola/r-index/build/ri-locate"
#------------------------Data folder--------------------------------------
filedirname = "/home/davide/data/data_r_index/"
datadirname = "/home/davide/data/"

# Logs
logdir_path = os.path.join(dirname,"logs/")

testRINDEX_log_path = os.path.join(logdir_path,"testRINDEX/" + datetime.datetime.now().strftime("%Y%m%d"))

logdirs = [testRINDEX_log_path] 

timeout = 86400 * 2 # 24h * 2

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    print("Creating the log directories")

    # create the logdir if not exists
    for logdir in logdirs:
        if not os.path.exists(logdir):
            os.makedirs(logdir)

    # get main bigbwt directory
    args.testdir = os.path.split(sys.argv[0])[0]
    # parameters
    datasets = [["salmonella_assemblies.fasta","salmonella_plain.fasta",32]]
    len_patt = [100,1000,10000]
    no_patt =  [100000,250000,400000,550000,700000,850000,1000000]

    running = {
        "er-index-count": True,
        "er-index-locate": True,
        "r-indexNicola-count": True,
        "r-indexNicola-locate": True,
        "r-indexTaher": False
    }

    print("Check the r-indexes construction and query times.")

    csvdir_path = os.path.join(logdir_path,"csvs/")
    if not os.path.exists(csvdir_path):
        os.makedirs(csvdir_path)

    csv_filename = os.path.join(csvdir_path, datetime.datetime.now().strftime("%Y%m%d") + ".csv")

    today = datetime.datetime.now().strftime("%Y%m%d") 

    print("Sending summary to: " + csv_filename)

    with open(csv_filename,"w") as csv_file:
        # Write CSV header 
        csv_file.write("Datastructure,Nseq,SeqLen,Dataset,TotLength,Query,PattLen,NoPatt,Time,MemPeak,CPU\n")

        # ---------- parsing of the input file
        start0 = start1 = start = time.time()
        for x in range(len(datasets)):
            for j in range(len(len_patt)):
                for y in range(len(no_patt)):

                    dname = datasets[x][0]
                    dname2 = datasets[x][1]
                    plen = len_patt[j]
                    nopat = no_patt[y]
                    
                    print("Computing count and locate queries of: " + dname + " and " + dname2)
                    print("Pattern length: " + str(plen) + "\nNumber of patterns: " + str(nopat))

                    start1 = time.time()

                    filename = os.path.join(filedirname,dname)
                    filename2 = os.path.join(filedirname,dname2)

                    # computing file stats using seqkit
                    print("Processing file " + filename)
                    command = "seqkit stats " + filename + " | tail -n 1"
                    res = subprocess.check_output(command, shell=True)
                    dataname = list(filter(('').__ne__, res.decode().split(" ")))[0]
                    nseq = int(list(filter(('').__ne__, res.decode().split(" ")))[3].replace(",", ""))
                    totlen = int(list(filter(('').__ne__, res.decode().split(" ")))[4].replace(",", ""))
                    avglen = float(list(filter(('').__ne__, res.decode().split(" ")))[6].replace(",", ""))
                    totlen2 = totlen + nseq
                    
                    if running["er-index-count"]:
                        # create pattern file
                        print("Computing pattern file of " + filename)
                        pattern_file = filename+"_"+str(plen)+"_"+str(nopat)+".pat"
                        command = "{exe} {input} {len} {nop} {out} 1 0".format(exe=genpattern_exe, input=filename, len=plen, nop=nopat,
                                                                               out=pattern_file)
                        print("=== " + command)
                        subprocess.check_output(command, shell=True)
                        # create count command 
                        if( datasets[x][2] == 64 ):
                            command = "{exe} {input} -q 0 -b 2 -f -v -p {pfile}".format(exe=er_index_exe_64, input=filename, pfile=pattern_file)
                            print(command)
                            csv_preamble = "eBWT r-index," + str(nseq) + "," + str(avglen) + "," + filename + "," + str(totlen) + "," + "count" + "," + str(plen) + "," + str(nopat) + ","
                        else:
                            command = "{exe} {input} -q 0 -b 2 -f -v -p {pfile}".format(exe=er_index_exe, input=filename, pfile=pattern_file)
                            print(command)
                            csv_preamble = "eBWT r-index," + str(nseq) + "," + str(avglen) + "," + filename + "," + str(totlen) + "," + "count" + "," + str(plen) + "," + str(nopat) + ","

                        running["er-index-count"] = test(command, dname, filename, testRINDEX_log_path, csv_file, csv_preamble)

                    if running["r-indexNicola-count"]:
                        # create pattern file
                        print("Computing pattern file of " + filename)
                        pattern_file = filename+"_"+str(plen)+"_"+str(nopat)+".pat"
                        command = "{exe} {input} {len} {nop} {out} 1 0".format(exe=genpattern_exe, input=filename, len=plen, nop=nopat,
                                                                               out=pattern_file)
                        print("=== " + command)
                        subprocess.check_output(command, shell=True)
                        # create count command 
                        command = "{exe} {index} {pfile}".format(exe=r_index_nic_count_exe, index=filename+".ri", pfile=pattern_file)
                        print(command)
                        csv_preamble = "r-index," + str(nseq) + "," + str(avglen) + "," + filename + "," + str(totlen) + "," + "count" + "," + str(plen) + "," + str(nopat) + ","

                        running["r-indexNicola-count"] = test(command, dname, filename, testRINDEX_log_path, csv_file, csv_preamble)

                    if running["er-index-locate"]:
                        # create pattern file
                        print("Computing pattern file of " + filename)
                        pattern_file = filename+"_"+str(plen)+"_"+str(nopat)+".pat"
                        command = "{exe} {input} {len} {nop} {out} 1 0".format(exe=genpattern_exe, input=filename, len=plen, nop=nopat,
                                                                               out=pattern_file)
                        print("=== " + command)
                        subprocess.check_output(command, shell=True)
                        # create count command 
                        if( datasets[x][2] == 64 ):
                            command = "{exe} {input} -q 1 -b 2 -f -v -p {pfile}".format(exe=er_index_exe_64, input=filename, pfile=pattern_file)
                            print(command)
                            csv_preamble = "eBWT r-index," + str(nseq) + "," + str(avglen) + "," + filename + "," + str(totlen) + "," + "count" + ","
                        else:
                            command = "{exe} {input} -q 1 -b 2 -f -v -p {pfile}".format(exe=er_index_exe, input=filename, pfile=pattern_file)
                            print(command)
                            csv_preamble = "eBWT r-index," + str(nseq) + "," + str(avglen) + "," + filename + "," + str(totlen) + "," + "locate" + "," + str(plen) + "," + str(nopat) + ","

                        running["er-index-locate"] = test(command, dname, filename, testRINDEX_log_path, csv_file, csv_preamble)

                    if running["r-indexNicola-locate"]:
                        # create pattern file
                        print("Computing pattern file of " + filename)
                        pattern_file = filename+"_"+str(plen)+"_"+str(nopat)+".pat"
                        command = "{exe} {input} {len} {nop} {out} 1 0".format(exe=genpattern_exe, input=filename, len=plen, nop=nopat,
                                                                               out=pattern_file)
                        print("=== " + command)
                        subprocess.check_output(command, shell=True)
                        # create count command 
                        command = "{exe} {index} {pfile}".format(exe=r_index_nic_locate_exe, index=filename+".ri", pfile=pattern_file)
                        print(command)
                        csv_preamble = "r-index," + str(nseq) + "," + str(avglen) + "," + filename + "," + str(totlen) + "," + "locate" + "," + str(plen) + "," + str(nopat) + ","

                        running["r-indexNicola-locate"] = test(command, dname, filename, testRINDEX_log_path, csv_file, csv_preamble)

                    print("Total length time: {0:.4f}".format(time.time()-start1), flush=True)

        print("Total experiment time: {0:.4f}".format(time.time()-start0), flush=True)
        print("==== Done")

def test(command, base_filename, filename, log_path, csv_file, csv_preamble):

    logfile_name = os.path.join(log_path, base_filename + ".log")
    print("\nSending logging messages to file:", logfile_name, flush=True)

    resfile_name = os.path.join(log_path, base_filename + ".res")
    print("Sending results messages to file:", resfile_name, flush=True)

    elapsed_time = 0.0

    # Clear page cache
    with open(filename,"r") as fd:
        length = os.path.getsize(filename)
        os.posix_fadvise(fd.fileno(), 0, length, os.POSIX_FADV_DONTNEED)

    with open(logfile_name,"a") as logfile:

        command = "/usr/bin/time -v -o {info} ".format(info=resfile_name) + command

        print("Command:", command, flush=True)
        start = time.time()
        ret = execute_command(command,logfile,logfile_name,timeout)
        elapsed_time = time.time()-start
        elapsed_time2 = 0
        print("Elapsed time: {0:.4f}".format(elapsed_time), flush=True)

        # ---- print elapsed time to file
        command = "echo Total construction time: {0:.4f}  \n\n".format(elapsed_time, flush=True)
        if(execute_command(command,logfile,logfile_name,)!=True):
            return

        rss = 0; cpu = 0;
        if ret:
            # Get resident set size from res file
            rss_o = subprocess.check_output("awk -F ':' '{{print $2}}' {log} | sed -n '10p'".format(log=resfile_name),shell=True)
            rss = int(re.search(r'\d+',rss_o.decode("utf-8")).group())
            # Get CPU time
            cpu_o = subprocess.check_output("awk -F ':' '{{print $2}}' {log} | sed -n '2p'".format(log=resfile_name),shell=True)
            cpu = cpu_o.decode("utf-8")
            cpu = cpu[1:-1]

        # complete csv line
        csv_preamble += str(elapsed_time) + "," + str(rss) + "," + str(cpu) + "\n"
        csv_file.write(csv_preamble)
        print(csv_preamble)

    return ret

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name,timeout=None,env=None):
    try:
        p = subprocess.Popen(command.split(),stdout=logfile,stderr=logfile,env=env,shell=False)
        # p = subprocess.Popen(command,stdout=logfile,stderr=logfile,env=env,shell=True)

        if timeout != None:

            try:
                p.communicate(timeout=timeout)
            except subprocess.TimeoutExpired as e:
                print("Timeout executing command line:")
                print("\t"+ command)
                # get all descendant pids
                pids_o = subprocess.check_output("pstree -p {pid} | grep -o '([0-9]\+)' | grep -o '[0-9]\+'".format(pid=p.pid),shell=True)
                pids = re.findall('\d+', pids_o.decode("utf-8"))
                print("Killing PIDs: " + "".join(str(pid) + " " for pid in pids))
                print("Check log file: " + logfile_name)
                for pid in pids:
                    os.kill(int(pid),signal.SIGKILL)
                    os.kill(int(pid),signal.SIGTERM)
                return False

    except subprocess.CalledProcessError as e:
        print("Error executing command line:")
        print("\t"+ command)
        print("Check log file: " + logfile_name)
        return False
    

    return True


if __name__ == '__main__':
    main()
