#!/usr/bin/env python
import os
from subprocess import check_call
import sys
import os.path

sys.path.append(os.path.abspath(os.path.curdir))
from shared import *

try:
    from params import *
except ImportError:
    print "params.py is not available"
    print "Try params_base.py"
    try:
        from params_base import *

    except ImportError:
        print "Try params_vel.py"
        from params_vel import *

gen_model_cords_bin = "/nesi/projects/nesi00213/tools/gen_model_cords"
#variables that cannot be found in params_base.py
CENTER_ORIGIN = '1'
GEOPROJ = '1'
DOCOORDS = '1'
BINDIR = gen_model_cords_bin #Although BINDIR was Null in original cshell, reassigned just in case BINDIR is used in some other scripts.

def main(outdir=None, debug=False):

    #assign variables according to cord_root    
    XLEN = float(nx)*float(hh)
    YLEN = float(ny)*float(hh)
    ZLEN = float(nz)*float(hh)

    if outdir is not None: #if this program is executed from genDomain, outdir is where all other vel model files will be placed.
        GRIDFILE = os.path.join(outdir, 'gridfile%s'%sufx)
        GRIDOUT = os.path.join(outdir, 'gridout%s'%sufx)
        MODEL_COORDS = os.path.join(outdir, 'model_coords%s'%sufx)
        MODELPARAMS = os.path.join(outdir, 'model_params%s'%sufx)
        MODEL_BOUNDS = os.path.join(outdir, 'model_bounds%s'%sufx)


    #generate Gridfile
    try:
        file_gridfile = open(GRIDFILE, "w")
    except IOError:
        print "cannot open file %s, please check file name or directory"%GRIDFILE
    else:
        file_gridfile.write("xlen=%f\n"%XLEN)
        file_gridfile.write("%10.4f %10.4f %13.6e\n"%(0.0,XLEN,float(hh)))
        file_gridfile.write("ylen=%f\n"%YLEN)
        file_gridfile.write("%10.4f %10.4f %13.6e\n"%(0.0,YLEN,float(hh)))
        file_gridfile.write("zlen=%f\n"%ZLEN)
        file_gridfile.write("%10.4f %10.4f %13.6e\n"%(0.0,ZLEN,float(hh)))
        file_gridfile.close()

        #run gen_model_cords by calling modul subprocess.call()

        if debug:
            errRedirect = ''
        else:
            errRedirect = '2>/dev/null'
        #calling executable "gen_model_cords"
        #using dictionary for better modification
        '''
        cmd = [gen_model_cords_bin, \
                'geoproj=' + GEOPROJ, 'gridfile=' + GRIDFILE, 'gridout=' + GRIDOUT, \
                'center_origin=' + CENTER_ORIGIN, 'do_coords=' + DOCOORDS, \
                'nzout=1', 'name=' + MODEL_COORDS, 'gzip=0', 'latfirst=0', \
                'modellon=' + MODEL_LON, 'modellat=' + MODEL_LAT, \
                'modelrot=' + MODEL_ROT, '1>'+ MODELPARAMS, errRedirect
                ]
        #exe(cmd,debug)
        print cmd
        check_call(cmd, shell=True)
        '''
        var_dic = {'BINDIR':BINDIR, 'GEOPROJ':GEOPROJ, 'GRIDFILE':GRIDFILE,'GRIDOUT':GRIDOUT,
                    'CENTER_ORIGIN':CENTER_ORIGIN,'DOCOORDS':DOCOORDS,
                    'MODEL_COORDS':MODEL_COORDS,
                    'MODEL_LON':MODEL_LON,'MODEL_LAT':MODEL_LAT,'MODEL_ROT':MODEL_ROT,'MODELPARAMS':MODELPARAMS,
                    'errRedirect': errRedirect
                   } 
        arg_string = ("{BINDIR} "
                   "geoproj={GEOPROJ} gridfile={GRIDFILE} gridout={GRIDOUT} "
                   "center_origin={CENTER_ORIGIN} do_coords={DOCOORDS} "
                   "nzout=1 name={MODEL_COORDS} gzip=0 latfirst=0 "
                   "modellon={MODEL_LON} modellat={MODEL_LAT} "
                   "modelrot={MODEL_ROT} 1> {MODELPARAMS} {errRedirect}").format(**var_dic)

        print arg_string
        check_call(arg_string,shell=True)
        
        """
        #legacy code for using call(). whole section has been remade into above
        #ver 1, long as code, impossible to modify
        #call("./gen_model_cords geoproj=%s gridfile=%s gridout=%s center_origin=%s do_coords=%s nzout=1 name=%s gzip=0 latfirst=0 modellon=%s modellat=%s modelrot=%s > %s"%(GEOPROJ,GRIDFILE,GRIDOUT,CENTER_ORIGIN,DOCOORDS,COORDFILE,MODEL_LON,MODEL_LAT,MODEL_ROT,PARAMFILE), shell = True)

        #ver2, call can not take multiple arguments "with multiple variables", hence these lines do not work
        check_call( ["./gen_model_cords gridfile=%s"%GRIDFILE, "geoproj=%s"%GEOPROJ ,"gridout=%s"%GRIDOUT,
               "center_origin=%s"%CENTER_ORIGIN, "do_coords=%s"%DOCOORDS, 
               "nzout=1", "name=%s"%COORDFILE, "gzip=0", "latfirst=0", 
               "modellon=%s"%MODEL_LON, "modellat=%s"%MODEL_LAT, 
               "modelrot=%s"%MODEL_ROT, ">", "%s"%PARAMFILE
              ], shell = True )
        #see codes above this section for modified version
        """
        #check
        if DOCOORDS =='1':
            #read MODEL_COORDS for read, and MODEL_BOUNDS for write
            try:
                file_coordfile = open(MODEL_COORDS, "r")
            except IOError:
                print "IO Eorror, cannot open files %s"%(MODEL_COORDS)
                sys.exit()
            try:
                file_boundfile = open(MODEL_BOUNDS, "w")
            except IOError:
                print "IO Eorror, cannot open files %s"%(MODEL_BOUNDS)
                sys.exit()

            #read lines from coordfiles
            while 1:
                line = file_coordfile.readline()
                if not line:
                    break
                if line[0] == '#':
                    continue
                if (float(line.split()[2]) == 0) or (float(line.split()[2]) == float(nx)-1) or (float(line.split()[3]) == 0) or (float(line.split()[3]) == float(ny)-1):
                    file_boundfile.write(line)
            file_coordfile.close()
            file_boundfile.close()
        #file_base.close()
        
if __name__ == "__main__":
    main()
