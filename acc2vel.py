from shared import *
import sys
import glob
import load_config as ldcfg
import os
qcore_cfg=ldcfg.load(os.path.join(os.path.dirname(os.path.realpath(__file__)),"qcore_config.json"))
int_bin = os.path.join(qcore_cfg['tools_dir'],'integ_diff') #int_bin = '/nesi/projects/nesi00213/tools/integ_diff'

if len(sys.argv) < 3:
    print "Usage:%s Acc_dir Vel_dir" %(sys.argv[0])
    sys.exit()

acc_path=os.path.relpath(sys.argv[1])
vel_path=os.path.abspath(sys.argv[2])

if not os.path.isdir(acc_path):
    print "No such Acc directory exists"
    sys.exit()
if not os.path.isdir(vel_path):
    os.makedirs(vel_path)
    print "Directory created: %s" %vel_path

files_in = glob.glob(os.path.join(acc_path,'*'))
files_out = [os.path.join(vel_path,os.path.basename(x)) for x in files_in]

for i,file_in in enumerate(files_in):
    file_out = files_out[i]
    print "file_in=%s \nfile_out=%s" %(file_in,file_out)

    exe([int_bin, 'integ=1', 'filein=%s' % file_in, 'inbin=0', 'outbin=0', 'fileout=%s' % file_out])

