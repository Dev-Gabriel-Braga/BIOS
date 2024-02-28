# Imports
from abaqus import *
from abaqusConstants import FULL, OFF
from odbAccess import *

# Initial Variables
base_name = 'ThreeBarTrussABAQUS'
output_precision = 8

# Submitting Job
job = mdb.JobFromInputFile(
    name = base_name,
    inputFileName = base_name + '.inp',
    nodalOutputPrecision = FULL    
)
job.submit(consistencyChecking=OFF)
job.waitForCompletion()

# Getting Stresses Values
odb = openOdb(base_name + '.odb')
stresses = odb.steps['Step-1'].frames[1].fieldOutputs['S'].values

# Writing Stresses in a Plain Text File
txt_file = open(base_name + '.txt', 'w')
format_string = "%." + str(output_precision - 1) + "e\n"
for s in stresses:
    txt_file.write(format_string % s.data[0])
txt_file.close()