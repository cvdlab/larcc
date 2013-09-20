from smplxn import *

#V,CV = simplexGrid([4,4,4])
V,CV = simplexGrid([1,1,1])
FV = facets(CV)
EV = facets(FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))

csrSignedBoundaryMat = signedBoundary (V,CV,FV)
print "\ncsrSignedBoundaryMat =\n", csrToMatrixRepresentation(csrSignedBoundaryMat)

boundaryCells_2 = signedBoundaryCells(V,CV,FV)
print "\nboundaryCells_2 =\n", boundaryCells_2

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
boundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells_2]

bndry = (V,boundaryFV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))
