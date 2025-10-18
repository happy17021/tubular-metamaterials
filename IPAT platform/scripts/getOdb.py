
#-*- coding:utf-8 -*-
import numpy as np
import sys
from odbAccess import *
import os
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import *
import displayGroupMdbToolset as dgm
import displayGroupOdbToolset as dgo
#########
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=199.322906494141,
    height=88.2349548339844)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
def mkdir(dir= "."):
	if(not os.path.exists(dir)):
		os.mkdir(dir)

mkdir("result")

odbName = sys.argv[-1]
outFrameCount = 10
# odbName = 'Job-1'
odb = session.openOdb(name=odbName + ".odb")
rgs = odb.steps['Step-1'].historyRegions
F_data = odb.steps['Step-1'].historyRegions[rgs.keys()[-1]].historyOutputs['RF2'].data
U_data = odb.steps['Step-1'].historyRegions[rgs.keys()[-1]].historyOutputs['U2'].data

data = []
for f,u in zip(F_data,U_data):
	data.append([abs(u[1]),abs(f[1]),0,0,0,0])
# out = sys.argv[-1] + ".txt"

np.savetxt(os.path.join("result",odbName + ".txt"),data,fmt='%f')
####
imgDir = os.path.join("result",odbName)
mkdir(imgDir)

def showVariable(variable = "SMISES"):
	try:
		if(variable.upper() == "SMISES"):
			session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
				variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT,
																				 'Mises'), )
	except BaseException as e:
		print(e)
		pass

#
def out2image(image_name = ""):
	print(image_name)
	# session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
	myViewport = session.viewports['Viewport: 1']
	# session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
	# myViewport.setValues(displayedObject=odb, width=190, height=890)
	session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
	# session.graphicsOptions.setValues(backgroudStyle=SOLID,backgroudColor='#FFFFFF')
	# session.printOptions.setValues(vpDecorations=OFF, reduceColors=False)
	session.printToFile(fileName=image_name, format=PNG, canvasObjects=(myViewport,))
	return

#

def changeFrame(stepId = 0,frameId =1):
	try:
		session.viewports['Viewport: 1'].odbDisplay.setFrame(step=stepId, frame=frameId)
	except BaseException as e:
		pass

####
steps = odb.steps
step = steps[steps.keys()[-1]]
stepId = len(steps)-1
frames = step.frames
#
N = len(frames)
#
showVariable("SMISES")
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    renderBeamProfiles=ON, renderShellThickness=ON)
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)
#############################
if(N <= outFrameCount):
	for i in range(N):
		changeFrame(stepId,i+1)
		session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
		out = os.path.join("result","%s//%s.png"%(odbName,i+1))
		out2image(out)
		# break
else:
	arr = np.linspace(0, N,outFrameCount)
	int_arr = np.round(arr).astype(int)
	for i,frameId in enumerate(int_arr):
		changeFrame(stepId, frameId)
		out = os.path.join("result", "%s//%s.png" % (odbName, i + 1))
		out2image(out)
		# break






