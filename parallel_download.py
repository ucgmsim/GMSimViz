hostname = 'fitzroy.nesi.org.nz' # remote hostname where SSH server is running
port = 22

glob_pattern='*.*'

import os
import glob
import paramiko
from mpi4py import MPI
import collections
import parallel_executor
import os.path
import humanize

import sys
import textwrap
import subprocess
import parallel_upload
import time
import stat

#get the size parameters from parallel_upload
terminal_cols = parallel_upload.terminal_cols
percentage_cols = parallel_upload.percentage_cols
progress_cols = parallel_upload.progress_cols
max_name_size = parallel_upload.max_name_size
humanized_size = parallel_upload.humanized_size
######
last_print = -1

def update_progress(msg,tag,msg_src):
    global last_print
    skip =False
    if last_print < 0:
        last_print = time.clock()
    now = time.clock()
    if int((now-last_print)*1000) < 10:
        skip = True
    else:
        last_print = now
        
    msg_src = msg_src -1
    filename, transfered, totalsize = msg
    fname_truncted = filename[:15]+'.'*3+filename[-42:]
    #if len(filename) > max_name_size:
    #    if len(os.path.basename(filename)) > max_name_size:
    #        filename = "..." + filename[(len(filename)-max_name_size):]
    #    else:
    #        filename = os.path.basename(filename)
    if (int(transfered) != int(totalsize)):
        output_msg = " "*(terminal_cols - percentage_cols)+"%.2f"%(float(transfered)*100.0 / float(totalsize))+'%  '
        output_msg = output_msg+ "\r"+" "*((terminal_cols - percentage_cols)-humanized_size)+"/"+humanize.naturalsize(totalsize,binary=True)
        output_msg = output_msg+ "\r"+" "*(((terminal_cols - percentage_cols)-humanized_size)-len(humanize.naturalsize(transfered,binary=True)))+humanize.naturalsize(transfered,binary=True)
        output_msg = output_msg+ "\r"+"Node %s: %s"%(str(msg_src+1),fname_truncted)
        output_dict[msg_src]['MSG'] = output_msg

    # revert cursor to the start
    if not skip:
        for i in output_dict:
            sys.stdout.write("\033[F")
            # extra revert because we printed a blank between progress and history
        sys.stdout.write("\033[F")

    # transfer completed
    if ( transfered == totalsize):
        if skip:
            for i in output_dict:
                sys.stdout.write("\033[F")
            sys.stdout.write("\033[F")
        completed = " Completed"
        # print from right to left, use \r to revert to start after every value printed
        print " "*(terminal_cols- len(completed)),completed,
        # print total size of file trasnfered
        print "\r"," "*((terminal_cols-len(completed))-len(humanize.naturalsize(totalsize,binary=True))),humanize.naturalsize(totalsize,binary=True),
        # print the file name
        print "\r", fname_truncted
        
        # update Node's status
        output_msg = " "*(terminal_cols - percentage_cols)+str(float(transfered)*100.0 / float(totalsize))+'%  '
        output_msg = output_msg+ "\r"+" "*((terminal_cols - percentage_cols)-humanized_size)+"/"+humanize.naturalsize(totalsize,binary=True)
        output_msg = output_msg+ "\r"+" "*(((terminal_cols - percentage_cols)-humanized_size)-len(humanize.naturalsize(transfered,binary=True)))+humanize.naturalsize(transfered,binary=True)
        output_msg = output_msg+ "\r"+"Node %s: %s"%(str(msg_src+1),fname_truncted)
        output_dict[msg_src]['MSG'] = output_msg
        if skip:
            print " "*terminal_cols
            for i in output_dict:
                print i['MSG']
                sys.stdout.flush()
    
    if not skip:
        # a blank line to split progress and history
        print " "*terminal_cols
        for i in output_dict:
            print i['MSG']
        sys.stdout.flush()

    

def download((local_file,remote_file)):
    def print_progress(bytes_transferred, bytes_total):
        comm.send([local_file, bytes_transferred,bytes_total],dest = 0, tag = 2)
    sftp.get(remote_file, local_file, print_progress)
    
def list_remote_files(sftp,directory):
    file_and_dir = sftp.listdir(directory)
    files = []
    dirs = [] 
    #print "start search in :",directory
    #print "found a list of :", file_and_dir 
    #determind if every element is a file or dir
    for x in file_and_dir:
        x_path = os.path.join(directory,x)
        #print "examine :",x_path
        if stat.S_ISDIR(sftp.stat(x_path).st_mode):
            #its dir
            dirs.append(x_path)
            #print "dir found: ",x_path
        else:
            files.append(x_path)
    
    for x in dirs:
        #print "recursive dir:", x
        new_files = list_remote_files(sftp,x)
        #print "find new files :", new_files
        files.extend(new_files)
    
    return files

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        print "SCP Download tool"

    dir_local, dir_remote, username_remote = parallel_upload.process_arguments(rank)
    #print "local ",dir_local
    #print "remote ", dir_remote
    
    sftp, t = parallel_upload.sftp_connect(hostname, port, username_remote)
    #if we are copying A/B/C to X/Y, we want A/B/C/* to be copied to X/Y/C/*
    if dir_remote.split(os.sep)[-1] == dir_local.split(os.sep)[-1]:
    #the destination is already X/Y/C, it's good. We just use the specified dir_remote
        pass
    else:
    #if destination is just X/Y, we are actually copying it to X/Y/C.
        dir_local=os.path.join(dir_local, dir_remote.split(os.sep)[-1])

    if rank == 0:
        try:
            os.makedirs(dir_local)
        except OSError, e:
            #print "%s already exists"%dir_local
            pass
    if rank == 0:
        print "Getting All files/directory info from Server, please hold..."
    files_remote = list_remote_files(sftp, dir_remote)
    dir_names = list(set([os.path.dirname(x) for x in files_remote]))
    dir_names = [os.path.relpath(x,dir_remote) for x in dir_names]
    files_local = [os.path.join(dir_local, os.path.relpath(x,dir_remote)) for x in files_remote]

    '''
    #print dir_remote
    if rank == 0:
        #print sftp.listdir(dir_remote)[0]
        #print stat.S_ISDIR(sftp.stat(os.path.join(dir_remote,sftp.listdir(dir_remote)[0])).st_mode)
        #files = list_remote_files(sftp, dir_remote)
        print files_remote

        print "local:",files_local
    '''
    if rank == 0:
        for x in dir_names:
            local_subdir = os.path.join(dir_local,x)
            try:
                os.makedirs(local_subdir)
            except OSError, e:
                #print "%s already exsits"%local_subdir
                pass
    #output_dict = []
        print "Node size :", size
        output_dict = []
        print "*"*40+"\n"
        for i in range(1, size):
            msg = "Node %s: Idle"%str(i)
            output_dict.append({"Node":i,"MSG":msg})
            print(output_dict[i-1]['MSG'])


    comm.Barrier()
    N=len(files_local)
    to_process=zip(files_local,files_remote)
    
    executor = parallel_executor.ParallelExecutor()
    executor.process_function_with_result(download,to_process,f_master = update_progress)

    t.close()





  
