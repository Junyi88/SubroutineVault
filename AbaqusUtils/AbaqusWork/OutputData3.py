# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

import sys
sys.path.insert(0,'Y://GitKrakenRes//AbaqusUtils//junyiPy')
sys.path.insert(0,'Y://GitKrakenRes//AbaqusUtils//abapy')

import numpy as np
from abapy.mesh import *
from abapy.postproc import *
from Rots import *

odb = session.openOdb(name= 'Y:\GitKrakenRes\AbaqusUtils\AbaqusWork\Exclude\MgO_10MicronC.odb' )
# odb = session.odbs['Y:\GitKrakenRes\AbaqusUtils\AbaqusWork\Exclude\MgO_10MicronC.odb']
frame = odb.steps[ 'Step-3' ].frames[-1]
dispField = frame.fieldOutputs['U'].values
my_part_instance = odb.rootAssembly.instances['PART-1-1']
numNodesTotal = len( my_part_instance.nodes )
# numNodesTotal = len( dispField )
# numElementsTotal = len( my_part_instance.elements )


label=[0]*numNodesTotal
x=[0]*numNodesTotal
y=[0]*numNodesTotal
z=[0]*numNodesTotal

for i in range( numNodesTotal ):
	curNode=my_part_instance.nodes[i]
	label[i]=curNode.label

	x[i]=curNode.coordinates[0]+dispField[i].data[0]
	y[i]=curNode.coordinates[1]+dispField[i].data[1]
	z[i]=curNode.coordinates[2]+dispField[i].data[2]
	# print i
nodes = Nodes(x=x,y=y,z=z, labels=label)
mesh = Mesh(nodes=nodes)

for i in range( numElementsTotal ):
	curElement =  list( my_part_instance.elements[i].connectivity )
	#print my_part_instance.elements[i].label
	mesh.add_element(label = my_part_instance.elements[i].label , 
		connectivity = curElement, space = 3 , name = 'hex8')
	print i



out = ''
out+=mesh.dump2vtk()
f = open("X2.vtk", "w")
f.write(out)
f.close()
# nodes = nodes(x=x,y=y,z=z, labels=label)
# mesh = mesh(nodes=nodes)

# out = ''
# out+=mesh.dump2vtk()
# f = open("x1.vtk", "w")
# f.write(out)
# f.close()