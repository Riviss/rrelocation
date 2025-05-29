#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 12:59:04 2025

@author: pgcseiscomp
"""

import pandas as pd

df = pd.read_csv('~/Documents/seismic_process/relocation_3D/nebc/data/aggregated/out.nllgrid3D.cat', delim_whitespace=True, header=None,
                 names=['year','month','day','hour','minute','second','event_id','latR','lonR','depR','mag','qID'
                        ,'cID','nbranch','qnpair','qndiffP','qndiffs','rmsP','rmsS','eh','ez','et','latC','lonC','depC'])

df.to_csv('~/Downloads/testy.csv')