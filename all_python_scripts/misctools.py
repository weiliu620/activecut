import sys
import os
import numpy as np

def (exp):
    """
    some basic test
    """

    exp = exp[:,1:]
    exp_log = log(exp)
    sample_mean = exp_log.mean(axis = 0)
    gene_mean = exp_log.mean(axis = 1)
