# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:07:41 2018

@author: oscar
"""

import matplotlib.pyplot as plt
import json

iters=[]
energy=[]

data = json.load(open('test.log'))

for iteration in data["Output"]:
    iters.append(iteration["Iteration"])
    energy.append(iteration["Energy"]["Mean"]+20.571809453)

plt.figure()        
plt.plot(iters,energy,color='red')
plt.title('Ground state Energy difference vs Iterations', bbox={'facecolor': '0.8', 'pad': 5})
plt.xlabel('Iteration')
plt.ylabel('Energy')
plt.savefig("ising2d_energydiff.png")
plt.show()

