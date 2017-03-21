#!/usr/bin/env python

## Copy files unattended over SSH using a glob pattern.
## It tries first to connect using a private key from a private key file
## or provided by an SSH agent. If RSA authentication fails, then 
## password login is attempted.

##
## DEPENDENT MODULES:
##      * paramiko, install it with `easy_install paramiko`

## NOTE: 1. The script assumes that the files on the source
##       computer are *always* newer that on the target;
##       2. no timestamps or file size comparisons are made
##       3. use at your own risk

hostname = 'fitzroy.nesi.org.nz' # remote hostname where SSH server is running
port = 22

glob_pattern='*.*'

import os
import glob
import paramiko
import md5
from mpi4py import MPI
import collections
import parallel_executor
import os.path
import humanize

import sys
#import getpass
#import curses
import textwrap
import subprocess
import time
#change this according to terminal size
terminal_cols = 100

#Only change these if the way of showing bitsize is changed
percentage_cols = 7
progress_cols = 20 + (percentage_cols+1)
max_name_size = terminal_cols - progress_cols -15
humanized_size = 11
####
last_print = -1

def update_progress(msg,tag,msg_src):
    global last_print
    skip = False
    if last_print < 0:
        last_print = time.clock()
    now = time.clock()
    if int((now-last_print)*1000) < 10:
        skip=True
    else:
        last_print = now

    # output_dicti[msg_src]
    msg_src = msg_src -1
    filename, transfered, totalsize = msg
    fname_truncted = filename[:15]+'.'*3+filename[-42:]
    #if len(filename) > max_name_size:
    #    if len(os.path.basename(filename)) > max_name_size:
    #        filename = "..." + filename[(len(filename)-max_name_size):]
    #    else:
    #        filename = os.path.basename(filename)

    if (int(transfered) != int(totalsize)):
        # update Node's status
        output_msg = " "*(terminal_cols - percentage_cols)+"%.2f"%(float(transfered)*100.0 / float(totalsize))+'%  '
        output_msg = output_msg+ "\r"+" "*((terminal_cols - percentage_cols)-humanized_size)+"/"+humanize.naturalsize(totalsize,binary=True)
        output_msg = output_msg+ "\r"+" "*(((terminal_cols - percentage_cols)-humanized_size)-len(humanize.naturalsize(    transfered,binary=True)))+humanize.naturalsize(transfered,binary=True)
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
        print "\r"," "*((terminal_cols-len(completed))-len(humanize.naturalsize(totalsize,binary=True)) ),humanize.naturalsize(totalsize,binary=True),
        # print the file name
        print "\r", fname_truncted
        
        # update Node's status
        output_msg = " "*(terminal_cols - percentage_cols)+str(float(transfered)*100.0 / float(totalsize))+'%  '
        output_msg = output_msg+ "\r"+" "*((terminal_cols - percentage_cols)-humanized_size)+"/"+humanize.naturalsize(totalsize,binary=True)
        output_msg = output_msg+ "\r"+" "*(((terminal_cols - percentage_cols)-humanized_size)-len(humanize.naturalsize(    transfered,binary=True)))+humanize.naturalsize(transfered,binary=True) 
        output_msg = output_msg+ "\r"+"Node %s: %s"%(str(msg_src+1),fname_truncted)
        output_dict[msg_src]['MSG'] = output_msg
        if skip:
            print " "*terminal_cols
            for i in output_dict:
                print i['MSG']
                # force flush stdout.
                sys.stdout.flush()
        
 
    #if not skip:
    # a blank line to split progress and history
    
    if not skip:
        print " "*terminal_cols

        # print out all Node's status
        for i in output_dict:
            print i['MSG']
        # force flush stdout.
        sys.stdout.flush() 
      

def agent_auth(transport, username,rsa_private_key):
    """
    Attempt to authenticate to the given transport using any of the private
    keys available from an SSH agent or from a local private RSA key file (assumes no pass phrase).
    """
#    print rsa_private_key
    try:
        ki = paramiko.RSAKey.from_private_key_file(rsa_private_key)
    except Exception, e:
        print 'Failed loading %s %s' % (rsa_private_key, e)

    agent = paramiko.Agent()
    agent_keys = agent.get_keys() + (ki,)
    if len(agent_keys) == 0:
        print "NO SSH keys found"
        return

    for key in agent_keys:
        #print 'Trying ssh-agent key %s' % key.get_fingerprint().encode('hex'),
        try:
            transport.auth_publickey(username, key)
            #print '... success!'
            return
        except paramiko.SSHException, e:
            print '... failed!', e
#            sys.exit()


def upload((local_file, remote_file)):
    def print_progress(bytes_transferred, bytes_total):
        #print '%-40s' % os.path.basename(local_file), humanize.naturalsize(bytes_transferred, binary=True), ' / ', humanize.naturalsize(bytes_total, binary=True), '%.2f' % (bytes_transferred * 100.0 / bytes_total), '%'
        #update progess by sending message(with tag=2) to master
        comm.send([local_file, bytes_transferred,bytes_total],dest = 0, tag = 2) 
#    local_file = os.path.join(dir_local, fname)
#    remote_file = dir_remote + '/' + os.path.basename(fname)

    #print 'Copying', local_file, 'to ', remote_file
    #print "file:",local_file
    #print "remote:",remote_file
    sftp.put(local_file, remote_file, print_progress)
    #print 'Completed: ', local_file, '--> ', remote_file
    return

def list_all_files(dir,files_only=[]):
    file_and_dir=glob.glob(dir+'/*')
    files = []
    dirs = []
    for x in file_and_dir:
        if os.path.isdir(x):
            dirs.append(x)
        else:
            files.append(x)
    for x in dirs:
        new_files = list_all_files(x)
        files.extend(new_files)

    return files




def process_arguments(rank):
    if len(sys.argv) != 4:
        if rank == 0:
            print "Usage: %s dir_local dir_remote username_remote" % sys.argv[0]
        sys.exit(0)
    dir_local = os.path.abspath(sys.argv[1])  # a training / will be removed
    dir_remote = os.path.abspath(sys.argv[2])
    # username_local = getpass.getuser()
    username_remote = sys.argv[3]
    return dir_local, dir_remote, username_remote


def sftp_connect(hostname, port, username_remote):
    rsa_private_key = os.path.join(os.getenv("HOME"), ".ssh/id_rsa")
    # print rsa_private_key
    # get host key, if we know one
    hostkeytype = None
    hostkey = None
    try:
        host_keys = paramiko.util.load_host_keys(os.path.expanduser('~/.ssh/known_hosts'))
    except IOError:
        try:
            # try ~/ssh/ too, e.g. on windows
            host_keys = paramiko.util.load_host_keys(os.path.expanduser('~/ssh/known_hosts'))
        except IOError:
            print '*** Unable to open host keys file'
            host_keys = {}

    # print host_keys.keys()
    if hostname in host_keys.keys():
        hostkeytype = host_keys[hostname].keys()[0]
        hostkey = host_keys[hostname][hostkeytype]
        # print 'Using host key of type %s' % hostkeytype

    # now, connect and use paramiko Transport to negotiate SSH2 across the connection
    # try:

    # dirlist on remote host
    #    dirlist = sftp.listdir('.')
    #    print "Dirlist:", dirlist
    print 'Establishing SSH connection to:', hostname, port, '...'
    t = paramiko.Transport((hostname, port))
    t.start_client()
    agent_auth(t, username_remote, rsa_private_key)
    if not t.is_authenticated():
        print 'RSA key auth failed! Set up your ssh key first'
        sys.exit(0)
    t.open_session()
    return paramiko.SFTPClient.from_transport(t), t

if __name__ == '__main__':
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    if rank == 0:
        print "SCP Upload tool"    

    dir_local, dir_remote, username_remote = process_arguments(rank)

    sftp, t = sftp_connect(hostname, port,username_remote)

    #if we are copying A/B/C to X/Y, we want A/B/C/* to be copied to X/Y/C/*
    if dir_local.split(os.sep)[-1] == dir_remote.split(os.sep)[-1]:
        #the destination is already X/Y/C, it's good. We just use the specified dir_remote
        pass
    else:
        #if destination is just X/Y, we are actually copying it to X/Y/C.
        dir_remote=os.path.join(dir_remote, dir_local.split(os.sep)[-1])
        #print dir_remotel

    if rank == 0:
        try:
            sftp.mkdir(dir_remote)
        except IOError, e:
            #print "%s already exists"%dir_remote
            pass

    files_local = list_all_files(dir_local)
    dir_names = list(set([os.path.dirname(x) for x in files_local]))
    dir_names = [os.path.relpath(x,dir_local) for x in dir_names] #relative paths
    files_remote = [os.path.join(dir_remote,os.path.relpath(x,dir_local)) for x in files_local]

    if rank == 0:
        #print "Local files:", files_local
        #print "Remote files:", files_remote
        for x in dir_names:
            remote_subdir = os.path.join(dir_remote,x)
            #split the directory and start mkdir from the root level to bottom
            split_dir = remote_subdir.split(os.sep)
            joined_dir = []
            if split_dir[0] == '':
                split_dir[0] = os.sep
            for i in  split_dir:
                joined_dir = os.path.join(joined_dir,i)
                try:
                    #print "mkdir %s" %remote_subdir
                    sftp.mkdir(joined_dir)
                except IOError, e:
                    #print "%s already exists"%joined_dir
                    pass
        print "Node size :",size
        output_dict = []
        #print output_dict
        print "*"*40+"\n"
        for i in range(1,size):
            msg = "Node %s: Idle"%str(i)
            output_dict.append({"Node":i,"MSG":msg})
            print(output_dict[i-1]['MSG'])
               

    comm.Barrier()
    N=len(files_local)
    to_process=zip(files_local,files_remote)

    executor = parallel_executor.ParallelExecutor()
    executor.process_function_with_result(upload,to_process,f_master = update_progress)
        
    t.close()

