# -*- coding: utf-8 -*-
"""
rmnoise.py

Function to filter out the noise region of my spectra in the different analysis notebook

@author: Aurelien Legoupil
"""

import numpy as np

def firstid(data, half_width_test = 40, noise_threshold = -75):
    extdata = np.zeros(len(data)+2*half_width_test)
    extdata[half_width_test:half_width_test+len(data)]=data
    TESTF_bool = [not(np.any(extdata[i:i+2*half_width_test] < noise_threshold)) for i in range(len(data))]
    try:
        I=next(id for id,bool in enumerate(TESTF_bool) if bool)
    except StopIteration:
        I=len(data)
    return I