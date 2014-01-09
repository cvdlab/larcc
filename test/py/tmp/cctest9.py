from simplexn import *
from larcc import *

V,CV = larSimplexGrid([5,5,5])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))

FV = larSimplexFacets(CV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))

EV = larSimplexFacets(FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))

np.set_printoptions(threshold='nan')

csrSignedBoundaryMat = signedBoundary (V,CV,FV)
Z = csrToMatrixRepresentation(csrSignedBoundaryMat)
print "\ncsrSignedBoundaryMat =\n", Z

from pylab import *
matshow(Z)
show()



boundaryCells_2 = signedBoundaryCells(V,CV,FV)
print "\nboundaryCells_2 =\n", boundaryCells_2

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
boundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells_2]

bndry = (V,boundaryFV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))
