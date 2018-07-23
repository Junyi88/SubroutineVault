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


label=[0]*numNodesTotal
x=[0]*numNodesTotal
y=[0]*numNodesTotal
z=[0]*numNodesTotal

f = open('XXX.txt','w')
f.write('PART-1-1  \n')
for i in range( numNodesTotal ):
	curNode=my_part_instance.nodes[i]
	f.write(str(i)+','+str(curNode.label)+'\n')

f.write('BERKOVICH-1 \n')

my_part_instance = odb.rootAssembly.instances['BERKOVICH-1']
numNodesTotal = len( my_part_instance.nodes )

for i in range( numNodesTotal ):
	curNode=my_part_instance.nodes[i]
	f.write(str(i)+','+str(curNode.label)+'\n')

#----------------------------------------------
numNodesTotal = len( frame.fieldOutputs['U'].values )	
f.write('x y z  \n')

f.close()

# numNodesTotal = len( dispField )
# numElementsTotal = len( my_part_instance.elements )


# session.odbs[name].steps[name].frames[i].fieldOutputs[name].values[i]\.instance.surfaces[name].nodes[i]


# nodes = nodes(x=x,y=y,z=z, labels=label)
# mesh = mesh(nodes=nodes)

# out = ''
# out+=mesh.dump2vtk()
# f = open("x1.vtk", "w")
# f.write(out)
# f.close()