# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

from abapy.mesh import *
from abapy.postproc import *

#outputDatabase = session.openOdb(name= 'Creep2.odb' )
odb = session.odbs['Y:/Tabassam1/Creep2.odb']
frame = odb.steps[ 'Step-1' ].frames[-1]
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

	x=x+[curNode.coordinates[0]]
	y=y+[curNode.coordinates[1]]
	z=z+[curNode.coordinates[2]]

nodeField = FieldOutput()
#print len(label)
#print len(x)
#print len(y)
#print len(z)

nodes = Nodes(x=x,y=y,z=z, labels=label)
mesh = Mesh(nodes=nodes)

for i in range( numElementsTotal ):
	curElement =  list( my_part_instance.elements[i].connectivity )
	#print my_part_instance.elements[i].label
	mesh.add_element(label = my_part_instance.elements[i].label , 
		connectivity = curElement, space = 3 , name = 'hex8')
#frame.fieldOutputs['U']
#frame.fieldOutputs['U'].values[0].data

u1 = FieldOutput()
for i in range( numNodesTotal ):
	curNode=my_part_instance.nodes[i].label
	outData = frame.fieldOutputs['U'].values[i].data[0]
	u1.add_data(data = outData, label = curNode)
	
U = VectorFieldOutput()
for i in range( numNodesTotal ):
	curNode=my_part_instance.nodes[i].label
	outData = frame.fieldOutputs['U'].values[i].data
	U.add_data(data1 = outData[0], data2 = outData[1], 
		data3 = outData[2], label = curNode)	
#elementField = frame.fieldOutputs['S']
#for i in range( numNodesTotal ):

PEEQ = FieldOutput(position='element')
for i in range( numElementsTotal  ):
	curElement =   my_part_instance.elements[i].label 
	outData = frame.fieldOutputs['S'].values[i].data[0]
	PEEQ.add_data(data = outData, label = curElement)

S = TensorFieldOutput(position='element')
for i in range( numElementsTotal  ):
	curElement =   my_part_instance.elements[i].label 
	outData = frame.fieldOutputs['S'].values[i].data
	S.add_data(data11 = outData[0], 
		data22 = outData[1], 
		data33 = outData[2], 
		data12 = outData[3], 
		data13 = outData[4], 
		data23 = outData[5], 
		label = curElement)
		
moron=GetTensorFieldOutput(odb, 'Step-1', -1, 'PART-1-1', 'element', 'S', labels=None, dti='I')
#print label
#print curElement
#print '===================================='

out = ''
out+=mesh.dump2vtk()
out+=u1.dump2vtk('u1')
out+=PEEQ.dump2vtk('PEEQ')
out+=U.dump2vtk('U')
out+=S.dump2vtk('S')
out+=moron.dump2vtk('moron')
f = open("FieldOutput-dump2vtk.vtk", "w")
f.write(out)
f.close()