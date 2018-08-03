# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:07:41 2018

@author: oscar
"""

import numpy as np
import matplotlib.pyplot as plt
import json

iters=[]
energy=[]

data = json.load(open('test.log'))

for iteration in data["Output"]:
    iters.append(iteration["Iteration"])
    energy.append(iteration["Energy"]["Mean"]+0.19910262)

plt.figure()        
plt.plot(iters[100:],energy[100:],color='red')
plt.title('Ground state Energy difference vs Iterations', bbox={'facecolor': '0.8', 'pad': 5})
plt.xlabel('Iteration')
plt.ylabel('Energy')
plt.savefig("Ryd2d_energydiff.png")
plt.show()
