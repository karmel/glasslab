'''
Created on Oct 15, 2012

@author: karmel
'''

def pandas_max(series, val=1):
    return series.apply(lambda x: max(x,val))

def pandas_min(series, val=1):
    return series.apply(lambda x: min(x,val))

def nonzero(series): return pandas_max(series, val=1)