from pyplasm import *
batches=[]
batches+=Batch.openObj("mergedMesh.obj")
octree=Octree(batches)
glcanvas=GLCanvas()
glcanvas.setOctree(octree)
glcanvas.runLoop()
