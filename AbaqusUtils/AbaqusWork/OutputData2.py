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

frame = odb.steps[ 'Step-3' ].frames[-1]
dispField = frame.fieldOutputs['U']
my_part_instance = odb.rootAssembly.instances['PART-1-1']
numNodesTotal = len( my_part_instance.nodes )
numElementsTotal = len( my_part_instance.elements )

label=[]
x=[]
y=[]
z=[]

for i in range( numNodesTotal ):
	curNode=my_part_instance.nodes[i]
	label=label+[curNode.label]

	x=x+[curNode.coordinates[0]+frame.fieldOutputs['U'].values[i].data[0]]
	y=y+[curNode.coordinates[1]+frame.fieldOutputs['U'].values[i].data[1]]
	z=z+[curNode.coordinates[2]+frame.fieldOutputs['U'].values[i].data[2]]

nodes = Nodes(x=x,y=y,z=z, labels=label)
mesh = Mesh(nodes=nodes)

out = ''
out+=mesh.dump2vtk()
f = open("X1.vtk", "w")
f.write(out)
f.close()