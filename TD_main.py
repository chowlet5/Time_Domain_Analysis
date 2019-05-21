import numpy as np
import pandas as pd
import yaml 

num_scheme = 'implicit'

input_file = 'TD_input.yml'

        
with open(input_file,'r') as stream:
    config = yaml.load(stream)['TD Input']

    n = config['num_floor']
    k = config['num_modes']
    T = config['duration']
    sf = config['sampling_freq']
N =T*sf
dt = 1/(sf/100)

