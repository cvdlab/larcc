from smplxn import *

V,FV = simplexGrid([10,10])
EV = simplexFacets(FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))

csrSignedBoundaryMat = signedBoundary (V,FV,EV)
print "\ncsrSignedBoundaryMat =\n", csrToMatrixRepresentation(csrSignedBoundaryMat)

boundaryCells_1 = signedBoundaryCells(V,FV,EV)
print "\nboundaryCells_1 =\n", boundaryCells_1

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
boundaryEV = [EV[-k] if k<0 else swap(EV[k]) for k in boundaryCells_1]

bndry = (V,boundaryEV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))
