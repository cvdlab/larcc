import scipy.misc, numpy
from numpy.random import randint
rows, columns = 100,100
rowSize, columnSize = 10,10

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
scipy.misc.imsave('./outfile.png', image_array)

noiseFraction = 0.1
noiseQuantity = rows*columns*noiseFraction
k = 0
while k < noiseQuantity:
   i,j = randint(rows),randint(columns)
   if image_array[i,j] == 0: image_array[i,j] = 255
   else: image_array[i,j] = 0
   k += 1
scipy.misc.imsave('./outfile.png', image_array)

from pyplasm import *
minPoint, maxPoint = (20,20), (40,30)
indexRanges = zip(minPoint,maxPoint)
window = CART([range(min,max) for min,max in indexRanges])

imageShape = [rows,columns]
d = len(imageShape)
weights = [PROD(imageShape[(k+1):]) for k in range(d-1)]+[1]
imageCochain = image_array.reshape(PROD(imageShape))
windowChain = [INNERPROD([index,weights]) for index in window]
segmentChain = [cell for cell in windowChain if imageCochain[cell]==255]

for cell in segmentChain: imageCochain[cell] = 127
image_array = imageCochain.reshape(imageShape)
scipy.misc.imsave('./outfile.png', image_array)

