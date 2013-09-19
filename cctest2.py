from smplxn import *
from scipy.spatial import Delaunay

#------------------------------------------------------------------
def signedBoundary (V,CV,FV):
	# compute the set of pairs of indices to [boundary face,incident coface]
	coo = boundary(CV,FV).tocoo()
	pairs = [[coo.row[k],coo.col[k]] for k,val in enumerate(coo.data) if val != 0]
	
	# compute the [face, coface] pair as vertex lists
	vertLists = [[FV[pair[0]], CV[pair[1]]]for pair in pairs]
	
	# compute two n-cells to compare for sign
	cellPairs = [ [list(set(coface).difference(face))+face,coface] for face,coface in vertLists]
	
	# compute the local indices of missing boundary cofaces
	missingVertIndices = [ coface.index(list(set(coface).difference(face))[0]) for face,coface in vertLists]
	
	# compute the point matrices to compare for sign
	pointArrays = [ [[V[k]+[1.0] for k in facetCell], [V[k]+[1.0] for k in cofaceCell]] for facetCell,cofaceCell in cellPairs]
	
	# signed incidence coefficients
	cofaceMats = TRANS(pointArrays)[1]
	cofaceSigns = AA(SIGN)(AA(np.linalg.det)(cofaceMats))
	faceSigns = AA(C(POWER)(-1))(missingVertIndices)
	signPairProd = AA(PROD)(TRANS([cofaceSigns,faceSigns]))
	
	# signed boundary matrix
	csrSignedBoundaryMat = csr_matrix( (signPairProd,TRANS(pairs)) )
	return csrSignedBoundaryMat

#------------------------------------------------------------------
def signedBoundaryCells(verts,cells,facets):
	csrBoundaryMat = signedBoundary(verts,cells,facets)
	csrTotalChain = totalChain(cells)
	csrBoundaryChain = matrixProduct(csrBoundaryMat, csrTotalChain)
	coo = csrBoundaryChain.tocoo()
	boundaryCells = list(coo.row * coo.data)
	return AA(int)(boundaryCells)

if __name__ == "__main__":
	
	V,CV = simplexGrid([4,4,4])
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

#------------------------------------------------------------------
verts = np.random.rand(1000, 3) # 1000 points in 3-d
verts = [AA(lambda x: 2*x)(VECTDIFF([vert,[0.5,0.5,0.5]])) for vert in verts]
verts = [vert for vert in verts if VECTNORM(vert) < 1.0]

tetra = Delaunay(verts)
cells = [cell for cell in tetra.vertices.tolist()
         if  ((verts[cell[0]][2]<0) and (verts[cell[1]][2]<0) and (verts[cell[2]][2]<0) and (verts[cell[3]][2]<0) ) ]

V, CV = verts, cells
FV = facets(CV)
#------------------------------------------------------------------
VIEW(MKPOL([V,AA(AA(lambda k:k+1))(FV),[]]))

boundaryCells_2 = boundaryCells(CV,FV)
print "\nboundaryCells_2 =\n", boundaryCells_2
bndry = (V,[FV[k] for k in boundaryCells_2])

VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))
#------------------------------------------------------------------


boundaryCells_2 = signedBoundaryCells(V,CV,FV)
print "\nboundaryCells_2 =\n", boundaryCells_2

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
boundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells_2]

bndry = (V,boundaryFV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))
