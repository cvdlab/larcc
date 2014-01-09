# -*- coding: utf-8 -*-
from simplexn import *
from largrid import *
from larcc import *

V,CV = larSimplexGrid([10,10,10])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))

FV = larSimplexFacets(CV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))

EV = larSimplexFacets(FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))

np.set_printoptions(threshold='nan')

csrSignedBoundaryMat = signedBoundary (V,CV,FV)
Z = csrToMatrixRepresentation(csrSignedBoundaryMat)
#print "\ncsrSignedBoundaryMat =\n", Z

from pylab import *
matshow(Z)
show()



boundaryCells_2 = signedBoundaryCells(V,CV,FV)
print "\nboundaryCells_2 =\n", boundaryCells_2

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
boundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells_2]

bndry = (V,boundaryFV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))


#---------------------------------------------------
boundaryEV = larSimplexFacets((V,boundaryFV), dim=2)[1]
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,boundaryEV))))

boundaryVertices = list(set(CAT(boundaryEV)))
print "boundaryVertices =",boundaryVertices,
bV = [V[k] for k in boundaryVertices]
print "boundaryV =",bV,

import collections
vertDict = collections.OrderedDict(zip(boundaryVertices,range(len(boundaryVertices))))
bEV = [[vertDict[e[0]], vertDict[e[1]]] for e in boundaryEV]
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
from pylab import *
matshow(Z)
show()


def Laplacian(G,W, x):
    Deltax = [[0,0,0]]*len(x)
    for i,j in (bEV + AA(REVERSE)(bEV)):
        Deltax[i] = VECTSUM([x[i], SCALARVECTPROD([W[i,j], VECTDIFF([x[i], x[j]]) ])])
    return Deltax

nV = AA(VECTSUM)(DISTR([bV,[-.5,-.5,-.5]]))

nV = Laplacian(bEV, W, nV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((nV, bEV))))
nV = Laplacian(bEV, W, nV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((nV, bEV))))
nV = Laplacian(bEV, W, nV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((nV, bEV))))




"""
    
Laplacian(G;W; x)
    new Δx = 0;
    for (e = (i,j) in E)
        Δxi = xi +w[i,j](xi - xj);
    end;
    return Δx;


LaplacianSmoothing(G;W;N;λ; x)
    new Δx
    for(i = 0 ; i < N ; i = i+1)
    Δx =Laplacian(G;W; x);
    x = x+λΔx;
    end;
    return;
    

TaubinSmoothing(G;W;N;λ;μ; x)
    new Δx
    for(i = 0 ; i < N ; i = i+1)
    Δx = Laplacian(G;W; x);
    if i is even
    x = x+λΔx;
    else
    x = x+μΔx;
    end;
    return;


FirFilter(G;W;N; f ; x)
    new x0
    = x
    new x1
    = Laplacian(G;W; x0
                );
    new x2
    = x0􀀀0:5x1
    new x = f0x0
    + f1x1
    for(i = 2 ; i < N ; i = i+1)
    x2
    =Laplacian(G;W; x1
               );
    x = x+ fix2;
    x0
    = x1;
    x1
    = x2;
    end;
    return;


"""