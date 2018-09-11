
import os

import json

config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.json")

with open(config_file, "r") as f:
    config = json.load(f)
