from smplxn import *

V = [[0.0, 0.0, 0.0],
	 [1.0, 0.0, 0.0],
	 [0.0, 1.0, 0.0],
	 [1.0, 1.0, 0.0],
	 [0.0, 0.0, 1.0],
	 [1.0, 0.0, 1.0],
	 [0.0, 1.0, 1.0],
	 [1.0, 1.0, 1.0]]

CV = [[0,1,2,3,4,5,6,7,8]]

FV = [[0,1,2,3],
	  [4,5,6,7],
	  [0,1,4,5],
	  [1,3,5,7],
	  [2,3,6,7],
	  [0,2,4,6]]

_,EV = larFacets((V,FV),dim=3)

VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))



csrSignedBoundaryMat = signedBoundary (V,CV,FV)
print "\ncsrSignedBoundaryMat =\n", csrToMatrixRepresentation(csrSignedBoundaryMat)

boundaryCells_2 = signedBoundaryCells(V,CV,FV)
print "\nboundaryCells_2 =\n", boundaryCells_2

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
boundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells_2]

bndry = (V,boundaryFV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))
