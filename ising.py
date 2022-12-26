# ising/ising.py : make input and output of 2d ising model

import os
import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pyising import draw

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--draw', type=int, default=None, help='<MCS> : DrawOutput')
args = parser.parse_args()

if args.draw:
	d = draw.Draw(args.draw)
	d.DrawOutput()
	sys.exit()
