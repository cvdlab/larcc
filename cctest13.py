from smplxn import *
import numpy as np

def LarIntervals(dom):
    def LarIntervals0(n):
        item = float(dom)/n
        ints = range(n+1)
        items = [item]*(n+1)
        verts = AA(LIST)(AA(PROD)(TRANS([ints,items])))
        cells = TRANS([ints[:-1],ints[1:]])
        return verts, cells
    return LarIntervals0

def boundCellsAdded(pairsOfBounds):
    def boundFilter0(model):
        verts,cells = model
        newcells = [[] for x in CAT(pairsOfBounds)]
        for k,v in enumerate(verts):
            print "k,v =",(k,v)
            for i,x in enumerate(v):
                print "i,x =",(i,x)
                bounds = pairsOfBounds[i]
                if x==bounds[0]: newcells[2*i] += [k]
                elif x==bounds[1]: newcells[2*i+1] += [k]
        return verts, cells+newcells
    return boundFilter0

mod_1 = LarIntervals(1)(4)
squares = larProduct([mod_1,mod_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
cubes = larProduct([squares,mod_1])
cubes = boundCellsAdded([[0,1],[0,1],[0,1]])(cubes)

V,CV = cubes
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV[:-6]))))
V,FV = larFacets((V,CV))
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
V,EV = larFacets((V,FV),dim=2)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))

CV = CV[:-6]
FV = [f for f in FV if len(f)==4]

csrBoundaryMat = boundary (CV,FV)
Z = csrToMatrixRepresentation(csrBoundaryMat)
#np.set_printoptions(threshold='nan')
#print "\ncsrSignedBoundaryMat =\n", Z
from pylab import *
matshow(Z)
show()

boundaryCells_2 = boundaryCells(CV,FV)
print "\nboundaryCells_2 =\n", boundaryCells_2

boundaryFV = [FV[k] for k in boundaryCells_2]
bndry = (V,boundaryFV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))


#---------------------------------------------------
boundaryEV = larFacets((V,boundaryFV), dim=2)[1]
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,boundaryEV))))

boundaryVertices = list(set(CAT(boundaryEV)))
print "boundaryVertices =",boundaryVertices,
bV = [V[k] for k in boundaryVertices]
print "boundaryV =",bV,

import collections
vertDict = collections.OrderedDict(zip(boundaryVertices,range(len(boundaryVertices))))
bEV = [[vertDict[e[0]], vertDict[e[1]]] for e in boundaryEV]
bFV = [[vertDict[v] for v in f] for f in boundaryFV]
bGraph = (bV, bEV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bGraph)))

csr_bEV = csrCreate(bEV)
Z = csrToMatrixRepresentation(csr_bEV)
from pylab import *
matshow(Z)
show()

VV = matrixProduct(csrTranspose(csr_bEV),csr_bEV).tocsr()
Z = csrToMatrixRepresentation(VV)
from pylab import *
matshow(Z)
show()

I,J,V = [],[],[]
D = VV.diagonal()
for k in range(VV.shape[0]):
    E = (VV[k]/D[k]).tocoo()
    V.extend(list(E.data)+[-1])
    J.extend(list(E.col)+[k])
    I.extend(len(E.data)*[k]+[k])
W = csr_matrix((V,(I,J)),VV.shape)
Z = csrToMatrixRepresentation(W)

np.set_printoptions(threshold='nan')
print "\ncsrSignedBoundaryMat =\n", Z

from pylab import *
matshow(Z)
show()


def Laplacian(G,W, x):
    Deltax = [[0,0,0]]*len(x)
    for i,j in G:
        Deltax[i] = VECTSUM([x[i], SCALARVECTPROD( [W[i,j], VECTDIFF([x[j], x[i]] ) ])])
    return Deltax

SCALARMATPROD = COMP([AA(SCALARVECTPROD),DISTL])
MATSUM = COMP([ AA(VECTSUM),TRANS ])

def LaplacianSmoothing(G,W, N,eta, x):
    for k in range(N):
        Deltax = Laplacian(G,W, x)
        x = MATSUM([ x, SCALARMATPROD([ eta, Deltax ]) ])
    return x


def TaubinSmoothing(G,W, N,eta,mu, x):
    for k in range(N):
        Deltax = Laplacian(G,W, x)
        if k % 2 == 0:
            x = MATSUM([ x, SCALARMATPROD([ eta, Deltax ]) ])
        else:
            x = MATSUM([ x, SCALARMATPROD([ mu, Deltax ]) ])
    return x


nV = AA(VECTSUM)(DISTR([bV,[-.5,-.5,-.5]]))

graph = bEV + AA(REVERSE)(bEV)

def mu(Lambda): return (3*Lambda -1)/(5*Lambda -3)

nV = LaplacianSmoothing(graph,W, 50,0.025, nV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((nV, bEV))))

nV = AA(VECTSUM)(DISTR([bV,[-.5,-.5,-.5]]))
nV = TaubinSmoothing(graph,W, 100,0.05,-0.03, nV)
VIEW(STRUCT(MKPOLS((nV, bFV))))


vmap = np.random.permutation(10)



"""
    
#Laplacian(G;W; x)
#    new Delta_x = 0;
#    for (e = (i,j) in E)
#        Delta_xi = xi +w[i,j](xi - xj);
#    end;
#    return Delta_x;
#
#
#LaplacianSmoothing(G;W;N;Lambda ; x)
#    new Delta_x
#    for(i = 0 ; i < N ; i = i+1)
#    Delta_x =Laplacian(G;W; x);
#    x = x+Lambda Delta_x;
#    end;
#    return;
#    
#
#TaubinSmoothing(G;W;N;Lambda ;Mu ; x)
#    new Delta_x
#    for(i = 0 ; i < N ; i = i+1)
#    Delta_x = Laplacian(G;W; x);
#    if i is even
#    x = x+Lambda Delta_x;
#    else
#    x = x+Mu Delta_x;
#    end;
#    return;
#
#
#FirFilter(G;W;N; f ; x)
#    new x0
#    = x
#    new x1
#    = Laplacian(G;W; x0
#                );
#    new x2
#    = x00:5x1
#    new x = f0x0
#    + f1x1
#    for(i = 2 ; i < N ; i = i+1)
#    x2
#    =Laplacian(G;W; x1
#               );
#    x = x+ fix2;
#    x0
#    = x1;
#    x1
#    = x2;
#    end;
#    return;


"""