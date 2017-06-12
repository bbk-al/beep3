#!/usr/bin/python
# -*- coding: utf-8 -*-
from pybeep import *

if __name__=="__main__":

    import sys
    from math import sqrt
    import constants

    cfg_filename = sys.argv[1]
    solver = BEEP(cfg_filename, False, True)


