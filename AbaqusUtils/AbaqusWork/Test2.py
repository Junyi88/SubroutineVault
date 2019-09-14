# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

from abapy.mesh import *
from abapy.postproc import *

odb = session.odbs['Y:/Tabassam1/Creep2.odb']
frame = odb.steps[ 'Step-1' ].frames[-1]
dispField = frame.fieldOutputs['U']
my_part_instance = odb.rootAssembly.instances['PART-1-1']
numNodesTotal = len( my_part_instance.nodes )
currNode=my_part_instance.nodes[1]


for i in range( numNodesTotal ):

	curNode =  my_part_instance.nodes[i]

	defNodePos = curNode.coordinates + dispField.values[i].data

	outFile.write( '\n' )

	for j in range( 3 ):

		outFile.write( str( defNodePos[j] ) + ' ' )

nodes = Nodes(x=x,y=y,z=z, labels=labels)
mesh = Mesh(nodes=nodes)
mesh.add_element(label = 1 , connectivity = [1,2,3], space = 2 , name = 'tri3') # triangle element
nodeField = FieldOutput()
nodeField.add_data(data = 0., label = 1)
nodeField.add_data(data = 10., label = 2)
nodeField.add_data(data = 20., label = 3)
elementField = FieldOutput(position='element')
elementField.add_data(label = 1, data =10.)
out = ''
out+=mesh.dump2vtk()
out+=nodeField.dump2vtk('nodeField')
out+=elementField.dump2vtk('elementField')
f = open("FieldOutput-dump2vtk.vtk", "w")
f.write(out)
f.close()