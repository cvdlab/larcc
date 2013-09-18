from smplxn import *

V,CV = simplexGrid([2,3,4])
grid_3d = (V,CV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_3d)))

SK2 = (V,facets(CV))
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK2)))

_,FV = SK2
SK1 = (V,facets(FV))
_,EV = SK1
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK1)))

