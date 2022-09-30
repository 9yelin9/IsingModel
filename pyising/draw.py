# ising/draw.py : draw diagrams of 2d ising model

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

class Draw:
	def __init__(self, MCS):
		self.MCS = MCS

		self.colors=['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
		self.markers = ['|', 'x', '*', 's']
		self.lss = ['solid', 'dashed', 'dashdot', 'dotted']
	
	def DrawPlot(self, data_type, ylim, ax):
		flist = [f for f in os.listdir('output/') if re.search('%s_MCS%s' % (data_type, self.MCS), f)]
		flist = sorted(flist, key=lambda x: int(re.sub('_L', '', re.search('_L\d+', x).group())))

		for i, fs in enumerate(flist):
			f = open('output/%s' % (fs), 'r')
			data = np.genfromtxt(f)
			f.close()

			if len(data):
				ax.plot(data[:, 0], data[:, 1], ls=self.lss[i], marker=self.markers[i], color=self.colors[i], label='L=%d' % int(re.sub('_L', '', re.search('_L\d+', fs).group())))

				ax.grid(True)
				ax.set_ylim(ylim)
				ax.set_xlabel('temperature')
				ax.set_ylabel('%s' % (data_type))
				ax.set_title('%s_MCS10e%s' % (data_type, self.MCS))
				ax.legend()

	def DrawOutput(self):
		fig, ax = plt.subplots(2, 2, figsize=(15, 8))

		self.DrawPlot("energy", [-2, -0.4], ax[0][0])
		self.DrawPlot("heat_capacity", [0, 1.6], ax[0][1])
		self.DrawPlot("abs_magnetization", [0, 1.05], ax[1][0])
		self.DrawPlot("mag_susceptibility", [0, 7], ax[1][1])

		fig.tight_layout()
		fig.savefig('ising.png')	
		plt.show()
