
import os
import sys

import json

config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
        'config.json')

with open(config_file, 'r') as f:
    qconfig = json.load(f)
