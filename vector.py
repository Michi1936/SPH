
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt


# In[41]:

import sys

readFile = open("analyzerTest.dat","r")
spaceCount=0
indexCount=0
pIndex=0
targetIndex=5
array =np.zeros((3,3))

for line in readFile:
    if line=='\n':
        spaceCount+=1
    if indexCount==targetIndex:
        print(line)
        splitted=line.split()
        for j in range(3):
            array[pIndex][j]=splitted[j]
        pIndex+=1
    if spaceCount==2:
        indexCount+=1
        spaceCount=0
    if indexCount>targetIndex:
        break
print(indexCount)


# In[42]:

array


# In[ ]:



