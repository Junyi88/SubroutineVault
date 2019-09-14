# create odb object from odb file
outputDatabase = session.openOdb(name= 'Creep2.odb' )

# get access to the nodal displacement data
frame = outputDatabase.steps[ 'Step-1' ].frames[-1]      

dispField = frame.fieldOutputs['U']

# get access to the part instance -- thru which u can access the undeformed nodal position coordinates
my_part_instance = outputDatabase.rootAssembly.instances['PART-1-1']


# Write deformed shape to vtk file 
# NOTE: if you want to export to a different file format, you should appropriate code for that file format here

outFile = open( 'deformed_shape.vtk' , 'w' )

    # write vtk header

outFile.write( '# vtk DataFile Version 3.0' )
outFile.write( '\nvtk output' )
outFile.write( '\nASCII' )
outFile.write( '\nDATASET UNSTRUCTURED_GRID' )

    # write points

numNodesTotal = len( my_part_instance.nodes )

outFile.write( '\n\nPOINTS ' + str( numNodesTotal ) + ' float' )

for i in range( numNodesTotal ):

	curNode =  my_part_instance.nodes[i]

	defNodePos = curNode.coordinates + dispField.values[i].data

	outFile.write( '\n' )

	for j in range( 3 ):

		outFile.write( str( defNodePos[j] ) + ' ' )

    # write cells

numElementsTotal = len( my_part_instance.elements )

outFile.write( '\n\nCELLS ' + str( numElementsTotal ) + ' ' + str( numElementsTotal * 8 ) )

for i in range( numElementsTotal ):

	curElement =  list( [4] + list( my_part_instance.elements[i].connectivity ) )

	outFile.write( '\n' )

	for j in range( 8 ):

		outFile.write( str( curElement[j] ) + ' ' )

    # write cell types

outFile.write( '\n\nCELL_TYPES ' + str( numElementsTotal ) )

for i in range( numElementsTotal ):

	outFile.write( '\n10' )

    # write cell data

outFile.write( '\n\nCELL_DATA ' + str( numElementsTotal ) )

    # write point data

outFile.write( '\n\nPOINT_DATA ' + str( numNodesTotal ) )

outFile.close()

outputDatabase.close()