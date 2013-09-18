# -*- coding: utf-8 -*-
"""
The MIT License
===============
    
Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from pyplasm import *
import collections
import scipy
import numpy as np
from scipy import zeros,arange,mat,amin,amax
from scipy.sparse import vstack,hstack,csr_matrix,coo_matrix,lil_matrix,triu

#------------------------------------------------------------------
#------------------------------------------------------------------
def CCOMB(vectors):
    """
        To create the convex combination of a list of vectors.
        Each vector is given as a list of coordinates.
        
        Return a vector.
        """
    return (COMP([ SCALARVECTPROD,CONS([ COMP([ DIV, CONS([K(1),LEN]) ]), VECTSUM ]) ]))(vectors)

def EXPLODE (sx,sy,sz):
    """
        To explode a HPC scene, given three real scaling parameters.
        sx,sy,sz >= 1.0
        
        Return a function to be applied to a list of HPC (Hierarchical Polyhedral Complex) objects.
        """
    def explode0 (scene):
        """
            To explode  a HPC scene, given as a list of HPC objects.
            Dimension-independent function (can be applied to points, edges, faces, cells, even mixed).
            Compute the centroid of each object, and apply to each of them a translation equal
            to the difference betwwen the scaled and the initial positions of its centroid.
            
            Return a single HPC object (the assembly of input objects, properly translated).
            """
        centers = [CCOMB(S1(UKPOL(obj))) for obj in scene]
        scalings = len(centers) * [S([1,2,3])([sx,sy,sz])]
        scaledCenters = [UK(APPLY(pair)) for pair in
                         zip(scalings, [MK(p) for p in centers])]
        translVectors = [ VECTDIFF((p,q)) for (p,q) in zip(scaledCenters, centers) ]
        translations = [ T([1,2,3])(v) for v in translVectors ]
        return STRUCT([ t(obj) for (t,obj) in zip(translations,scene) ])
    return explode0

def MKPOLS (model):
    """
        To MaKe a list of HPC objects from a LAR model.
        A LAR model is a pair, i.e. a Python tuple (V, FV), where
        -   V is the list of vertices, given as lists of coordinates;
        -   FV is the face-vertex relation, given as a list of faces,
        where each face is given as a list of vertex indices.
        
        Return a list of HPC objects.
        """
    V, FV = model
    pols = [MKPOL([[V[v] for v in f],[range(1,len(f)+1)], None]) for f in FV]
    return pols

#------------------------------------------------------------------
#------------------------------------------------------------------
def format(cmat,shape="csr"):
    """ Transform from list of triples (row,column,vale)
        to scipy.sparse corresponding formats.
        
        Return by default a csr format of a scipy sparse matrix.
        """
    n = len(cmat)
    data = arange(n)
    ij = arange(2*n).reshape(2,n)
    for k,item in enumerate(cmat):
        ij[0][k],ij[1][k],data[k] = item
    return scipy.sparse.coo_matrix((data, ij)).asformat(shape)

#------------------------------------------------------------------
def cooCreateFromBrc(ListOfListOfInt):
    COOm = [[k,col,1] for k,row in enumerate(ListOfListOfInt)
            for col in row ]
    return COOm

if __name__ == "__main__":
    print "\n>>> cooCreateFromBrc"
    V = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
    FV = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
    EV = [[0,1],[0,3],[1,2],[1,3],[1,4],[2,4],[2,5],[3,4],[4,5]]
    cooFV = cooCreateFromBrc(FV)
    cooEV = cooCreateFromBrc(EV)
    print "\ncooCreateFromBrc(FV) =\n", cooFV
    print "\ncooCreateFromBrc(EV) =\n", cooEV

#------------------------------------------------------------------
def csrCreateFromCoo(COOm):
    CSRm = format(COOm,"csr")
    return CSRm

if __name__ == "__main__":
    print "\n>>> csrCreateFromCoo"
    csrFV = csrCreateFromCoo(cooFV)
    csrEV = csrCreateFromCoo(cooEV)
    print "\ncsr(FV) =\n", repr(csrFV)
    print "\ncsr(EV) =\n", repr(csrEV)

#------------------------------------------------------------------
def csrCreate(BRCm,shape=(0,0)):
    if shape == (0,0):
        out = csrCreateFromCoo(cooCreateFromBrc(BRCm))
        return out
    else:
        CSRm = scipy.sparse.csr_matrix(shape)
        for i,j,v in cooCreateFromBrc(BRCm):
            CSRm[i,j] = v
        return CSRm

if __name__ == "__main__":
    print "\n>>> csrCreateFromCoo"
    V = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
    FV = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
    csrFV = csrCreate(FV)
    print "\ncsrCreate(FV) =\n", csrFV

#------------------------------------------------------------------
def csrGetNumberOfRows(CSRm):
    Int = CSRm.shape[0]
    return Int

if __name__ == "__main__":
    print "\n>>> csrGetNumberOfRows"
    print "\ncsrGetNumberOfRows(csrFV) =", csrGetNumberOfRows(csrFV)
    print "\ncsrGetNumberOfRows(csrEV) =", csrGetNumberOfRows(csrEV)

#------------------------------------------------------------------
def csrGetNumberOfColumns(CSRm):
    Int = CSRm.shape[1]
    return Int

if __name__ == "__main__":
    print "\n>>> csrGetNumberOfColumns"
    print "\ncsrGetNumberOfColumns(csrFV) =", csrGetNumberOfColumns(csrFV)
    print "\ncsrGetNumberOfColumns(csrEV) =", csrGetNumberOfColumns(csrEV)

#------------------------------------------------------------------
def csrToMatrixRepresentation(CSRm):
    nrows = csrGetNumberOfRows(CSRm)
    ncolumns = csrGetNumberOfColumns(CSRm)
    ScipyMat = zeros((nrows,ncolumns),int)
    C = CSRm.tocoo()
    for triple in zip(C.row,C.col,C.data):
        ScipyMat[triple[0],triple[1]] = triple[2]
    return ScipyMat

if __name__ == "__main__":
    print "\n>>> csrToMatrixRepresentation"
    print "\nFV =\n", csrToMatrixRepresentation(csrFV)
    print "\nEV =\n", csrToMatrixRepresentation(csrEV)

#------------------------------------------------------------------

def matrixProduct(CSRm1,CSRm2):
    CSRm = CSRm1 * CSRm2
    return CSRm

def csrTranspose(CSRm):
    CSRm = CSRm.T
    return CSRm

#------------------------------------------------------------------
def csrBoundaryFilter(CSRm, facetLengths):
    maxs = [max(CSRm[k].data) for k in range(CSRm.shape[0])]
    inputShape = CSRm.shape
    coo = CSRm.tocoo()
    for k in range(len(coo.data)):
        if coo.data[k]==maxs[coo.row[k]]: coo.data[k] = 1
        else: coo.data[k] = 0
    mtx = coo_matrix((coo.data, (coo.row, coo.col)), shape=inputShape)
    out = mtx.tocsr()
    return out

if __name__ == "__main__":
    print "\n>>> csrBoundaryFilter"
    csrEF = matrixProduct(csrFV, csrTranspose(csrEV)).T
    facetLengths = [csrCell.getnnz() for csrCell in csrEV]
    CSRm = csrBoundaryFilter(csrEF, facetLengths).T
    print "\ncsrMaxFilter(csrFE) =\n", csrToMatrixRepresentation(CSRm)

#------------------------------------------------------------------
# model = (vertices, topology)

def larProduct(models):
    model1,model2 = models
    V, cells1 = model1
    W, cells2 = model2
    verts = collections.OrderedDict(); k = 0
    for v in V:
        for w in W:
            vertex = tuple(v+w)
            if not verts.has_key(vertex):
                verts[vertex] = k
                k += 1
    cells = [ sorted([verts[tuple(V[v]+W[w])] for v in c1 for w in c2])
             for c1 in cells1 for c2 in cells2]
    
    model = AA(list)(verts.keys()), sorted(cells)
    return model

if __name__ == "__main__":
    geom_0,topol_0 = [[0.],[1.],[2.],[3.],[4.]],[[0,1],[1,2],[3,4]]
    geom_1,topol_1 = [[0.],[1.],[2.]], [[0,1],[1,2]]
    mod_0 = (geom_0,topol_0)
    mod_1 = (geom_1,topol_1)
    
    squares = larProduct([mod_0,mod_1])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
    
    cubes = INSL(larProduct)([mod_0,mod_1,mod_0]) # ==
    cubes = larProduct([squares,mod_0])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))

#------------------------------------------------------------------


if __name__ == "__main__":
    
    pass

