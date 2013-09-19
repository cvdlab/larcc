from smplxn import *

V,CV = simplexGrid([1,1,1])
#V,CV = simplexGrid([10,10,3])
grid_3d = (V,CV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_3d)))

SK2 = (V,facets(CV))
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK2)))

_,FV = SK2
SK1 = (V,facets(FV))
_,EV = SK1
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK1)))

csrC = totalChain(CV)
csrF = totalChain(FV)
csrE = totalChain(EV)

"""
print "\ntotalChain_0 =\n", csrToMatrixRepresentation(csrC)
print "\ntotalChain_1 =\n", csrToMatrixRepresentation(csrF)
print "\ntotalChain_2 =\n", csrToMatrixRepresentation(csrE)
"""

print "len(csrC), len(csrF), len(csrE) =",csrC.getnnz(),csrF.getnnz(),csrE.getnnz(),"\n"


def boundaryCells(cells,facets):
	csrBoundaryMat = boundary(cells,facets)
	csrChain = totalChain(cells)
	csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
	for k,value in enumerate(csrBoundaryChain.data):
		if value % 2 == 0: csrBoundaryChain.data[k] = 0
	boundaryCells = [k for k,val in enumerate(csrBoundaryChain.data.tolist()) if val == 1]
	return boundaryCells
	
	
boundaryCells_2 = boundaryCells(CV,FV)
#boundaryCells_1 = boundaryCells([FV[k] for k in boundaryCells_2],EV)

print "\nboundaryCells_2 =\n", boundaryCells_2
#print "\nboundaryCells_1 =\n", boundaryCells_1

boundary = (V,[FV[k] for k in boundaryCells_2])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(boundary)))

