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
numNodesTotal=len(frame.fieldOutputs['U'].values)

my_part_instance = odb.rootAssembly.instances['PART-1-1']
numNodesPart = len( my_part_instance.nodes )
numElementsPart = len( my_part_instance.elements )

# Construct Parts
label=[0]*numNodesPart
x=[0]*numNodesPart
y=[0]*numNodesPart
z=[0]*numNodesPart

#session.odbs[name].steps[name].frames[i].fieldOutputs[name].values[i]\.instance.surfaces[name].nodes[i]
# Run the loop
count=0
for i in range( numNodesTotal ):
	if (dispField [i].instance.name == 'PART-1-1'):
		label[count]=dispField [i].nodeLabel
		x[count]=dispField[i].instance.nodes[count].coordinates[0]+dispField[i].data[0]
		y[count]=dispField[i].instance.nodes[count].coordinates[1]+dispField[i].data[1]
		z[count]=dispField[i].instance.nodes[count].coordinates[2]+dispField[i].data[2]
		count=count+1
	print i

nodes = Nodes(x=x,y=y,z=z, labels=label)
mesh = Mesh(nodes=nodes)

# Now Add Values
SField = frame.fieldOutputs

count=0
SDV = TensorFieldOutput(position='element')
RR = VectorFieldOutput(position='element')

for i in range( len(SField['SDV1'].values)  ):
	curElement =  list( my_part_instance.elements[i].connectivity )
	#print my_part_instance.elements[i].label
	mesh.add_element(label = my_part_instance.elements[i].label , 
		connectivity = curElement, space = 3 , name = 'hex8')
	curElement =   SField['SDV199'].values[i].elementLabel
	SDV.add_data(
		data11 = SField['SDV199'].values[i].data, 
		data22 = SField['SDV203'].values[i].data, 
		data33 = SField['SDV207'].values[i].data, 
		data12 = SField['SDV200'].values[i].data, 
		data13 = SField['SDV201'].values[i].data, 
		data23 = SField['SDV204'].values[i].data, 
		label = curElement)
	data11 = SField['SDV199'].values[i].data
	data22 = SField['SDV203'].values[i].data
	data33 = SField['SDV207'].values[i].data
	data12 = SField['SDV200'].values[i].data
	data13 = SField['SDV201'].values[i].data
	data23 = SField['SDV204'].values[i].data
	R=np.array([[data11,data12,data13],[data12,data22,data23],
		[data13,data23,data33]])
	Ro=rotationMatrixToEulerAngles(R)
	RR.add_data(data1 = Ro[0], data2 = Ro[1], 
		data3 = Ro[2], label = curElement)
	print i

# Calculate Euler Angles

		
out = ''
out+=mesh.dump2vtk()
out+=SDV.dump2vtk('RotatinMat')
out+=RR.dump2vtk('Rot1')
f = open("X8.vtk", "w")
f.write(out)
f.close()