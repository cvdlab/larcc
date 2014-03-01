""" LAR implementation of morphological operators on multidimensional images."""

import scipy.misc, numpy
from numpy.random import randint
from pyplasm import *

""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')

import largrid
from largrid import *


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

def mapTupleToInt(shape):
   d = len(shape)
   weights = [PROD(shape[(k+1):]) for k in range(d-1)]+[1]
   
   def mapTupleToInt0(tuple):
      return INNERPROD([tuple,weights])
   return mapTupleToInt0

def setMaskWindow(window,image_array):
   minPoint, maxPoint = window
   imageShape = list(image_array.shape)
   indexRanges = zip(minPoint,maxPoint)
   tuples = CART([range(min,max) for min,max in indexRanges])
   
   imageCochain = image_array.reshape(PROD(imageShape))
   mapping = mapTupleToInt(imageShape)
   windowChain = [mapping(tuple) for tuple in tuples]
   segmentChain = [cell for cell in windowChain if imageCochain[cell]==255]
   
   for cell in segmentChain: imageCochain[cell] = 127
   image_array = imageCochain.reshape(imageShape)
   scipy.misc.imsave('./outfile.png', image_array)
   
   return segmentChain

def larImage(shape):
   """ Compute vertices and skeletons of an image of given shape """
   imageVerts,_ = larCuboids(list(shape))
   skeletons = gridSkeletons(list(shape))
   return imageVerts, skeletons

def boundaryOps(skeletons):
   """ CSR matrices of boundary operators from list of skeletons """
   return [boundary(skeletons[k+1],faces) 
      for k,faces in enumerate(skeletons[:-1])]

def visImageChain (shape,chain):
   imageVerts, skeletons = larImage(shape)
   chainLAR = [cell for k,cell in enumerate(skeletons[-1]) if k in chain]
   return imageVerts,chainLAR

def imageChainBoundary(shape):
   imageVerts, skeletons = larImage(shape)
   operators = boundaryOps(skeletons)
   cellNumber = PROD(list(shape))
   
   def imageChainBoundary0(k):
      csrBoundaryMat = operators[-1]
      facets = skeletons[k-1]
      
      def imageChainBoundary1(chain):
         csrChain = scipy.sparse.csr_matrix((cellNumber,1))
         for h in chain: csrChain[h,0] = 1
         csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
         for h,value in enumerate(csrBoundaryChain.data):
            if MOD([value,2]) == 0: csrBoundaryChain.data[h] = 0
         cooBoundaryChain = csrBoundaryChain.tocoo()
         boundaryCells = [cooBoundaryChain.row[h] 
            for h,val in enumerate(cooBoundaryChain.data) if val == 1]
         
         boundaryChainModel = imageVerts, [facets[h] for h in boundaryCells]     
         return boundaryChainModel
      
      return imageChainBoundary1
   return imageChainBoundary0

