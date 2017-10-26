
import sys
import os
import json

def load(config_file):
    #directory = os.path.dirname(os.path.abspath(__file__))
    #config_file = os.path.join(directory, "workflow_config.json")
    #print directory
    #exit(1)
    try:
        with open(config_file) as f:
            config_dict = json.load(f)
            return config_dict

    except IOError:
        print "cannot reach config file: %s" %config_file 
        print "This is a fatal error. Please contact someone from the software team."
        exit(1)



