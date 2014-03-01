
import scipy.misc, numpy
from numpy.random import randint
from pyplasm import *

""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')

import largrid
from largrid import *


import morph
from morph import *
 
rows, columns = 100,100
rowSize, columnSize = 10,10
shape = (rows, columns)
structure = (rowSize, columnSize)
image_array = randomImage(shape, structure, 0.3)
minPoint, maxPoint = (20,20), (40,30)
window = minPoint, maxPoint
segmentChain = setMaskWindow(window,image_array)
   
solid = visImageChain (shape,segmentChain)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(solid)))
b_rep = imageChainBoundary(shape)(2)(segmentChain)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(b_rep)))
