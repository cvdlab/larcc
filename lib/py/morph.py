""" LAR implementation of morphological operators on multidimensional images."""

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
 

def randomImage(shape, structure, noiseFraction=0.1):
   """ Generation of random image of given shape and structure. 
      Return a scipy.ndarray(shape)
   """
   rows, columns = shape
   rowSize, columnSize = structure                 
   random_array = randint(0, 255, size=(rowSize, columnSize))
   image_array = numpy.zeros((rows, columns))
   for i in range(rowSize):
      for j in range(columnSize):
         for h in range(i*rowSize,i*rowSize+rowSize): 
            for k in range(j*columnSize,j*columnSize+columnSize):
               if random_array[i,j] < 127:
                  image_array[h,k] = 0 
               else: 
                  image_array[h,k] = 255
   
   noiseQuantity = rows*columns*noiseFraction
   k = 0
   while k < noiseQuantity:
      i,j = randint(rows),randint(columns)
      if image_array[i,j] == 0: image_array[i,j] = 255
      else: image_array[i,j] = 0
      k += 1
   scipy.misc.imsave('./outfile.png', image_array)
   
   return image_array

def setMaskWindow(window,image_array):
   minPoint, maxPoint = window
   imageShape = list(image_array.shape)
   indexRanges = zip(minPoint,maxPoint)
   tuples = CART([range(min,max) for min,max in indexRanges])
   
   d = len(imageShape)
   weights = [PROD(imageShape[(k+1):]) for k in range(d-1)]+[1]
   imageCochain = image_array.reshape(PROD(imageShape))
   windowChain = [INNERPROD([index,weights]) for index in tuples]
   segmentChain = [cell for cell in windowChain if imageCochain[cell]==255]
   
   for cell in segmentChain: imageCochain[cell] = 127
   image_array = imageCochain.reshape(imageShape)
   scipy.misc.imsave('./outfile.png', image_array)
   
   return segmentChain

def visImageChain (shape,chain):
   imageShape = list(shape)
   model = larCuboids(imageShape)
   imageVerts = model[0]
   imageLAR = model[1]
   chainLAR = [cell for k,cell in enumerate(imageLAR) if k in chain]
   return imageVerts,chainLAR

def visImageChainBoundary (shape,chain):
   imageShape = list(shape)
   model = larCuboids(imageShape)
   imageVerts = model[0]
   skeletons = gridSkeletons(imageShape)
   facets = skeletons[-2]
   csrBoundaryMat = gridBoundaryMatrices(imageShape)[-1]
   csrChain = scipy.sparse.csr_matrix((PROD(imageShape),1))
   for k in chain: csrChain[k,0] = 1
   csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
   for k,value in enumerate(csrBoundaryChain.data):
      if MOD([value,2]) == 0: csrBoundaryChain.data[k] = 0
   cooBoundaryChain = csrBoundaryChain.tocoo()
   boundaryCells = [cooBoundaryChain.row[k] 
      for k,val in enumerate(cooBoundaryChain.data) if val == 1]
   return imageVerts,[facets[k] for k in boundaryCells]

