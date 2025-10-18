# -*- coding:utf-8 -*-
import sys
from abaqus import *
from abaqusConstants import *
import math
import numpy as np
# data preprocessing
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=158.338531494141,
    height=92.4791641235352)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
#####################################

def map_polyline_to_circle(points,R = 10.):
    # Step 3: map each point to circle
    total_length = R * math.pi * 2.
    circle_points = []
    for i,point in enumerate(points):
        theta = (point[0] / total_length) * 2. * math.pi
        x = R * math.cos(theta)
        y = R * math.sin(theta)
        circle_points.append((x,point[1],y))
    return circle_points
###############################################################
Mdb()
# Set user-defined parameters
grid = [3,3,3]
pix_width = 128.
pix_height = 128.
step_time = 0.5
meshSize = 1.0 
pattern = [6,3]
block = [12.8,12.8]
layers = 4
depth = 10.
radio = 1.*block[0]/pix_width
displacement = 50.
Thickness_cell = 1.5

NNNN = 200
# model creation
# N = 12


N = int(sys.argv[-1])
baseName = str(sys.argv[-2])

# baseName = "221"
# N = int(sys.argv[-1])
# baseName = str(sys.argv[-2])

# N = 4
# baseName = "RA_4"
N = int(sys.argv[-1])
baseName = str(sys.argv[-2])

block = np.asarray(block) #+ np.asarray([.15,.15])
L = grid[0] * block[0]*2.
W = grid[1] * block[0]*2.
H = pattern[1] * block[0]*2.

mat_data = {
    'e':36210,  
    'u':0.31,
    'density':2.65e-09,
    'plastic':[(0.1000, 0.00),
               (35.191, 0.01),
               (80.388, 0.02),
               (117.69, 0.03),
               (136.78, 0.04),
               (148.41, 0.05),
               (155.57, 0.06),
               (160.58, 0.07),
               (164.03, 0.08),
               (166.68, 0.09),
               (170.29, 0.10),
               (171.47, 0.11),
               (171.47, 0.12),
               (172.70, 0.13),
               (173.72, 0.14),
               (174.40, 0.15),
               (175.08, 0.16),
               (175.89, 0.17),
               (176.37, 0.18),
               (176.59, 0.19),
               (175.28, 0.20)]
}





layerD = 0.1
cutParts= []
for i in range(N):
    try:
        data = np.loadtxt('data/{}_{}.txt'.format(baseName,i+1))#*(radio+0.001)
        approx = data.tolist()
        approx.append(data[0])
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=200.0)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)

        # s.Spline(points=approx)
        for j in range(len(approx)-1):
            pt1,pt2 = approx[j],approx[j+1]
            s.Line(point1=pt1, point2=pt2)

        session.viewports['Viewport: 1'].view.setValues(width=135.349, height=74.4681,
            cameraPosition=(0.406673, 1.59753, 188.562), cameraTarget=(0.406673,
            1.59753, 0))
        p = mdb.models['Model-1'].Part(name='Part-{}'.format(i+1), dimensionality=THREE_D,
                                       type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts['Part-{}'.format(i+1)]
        cutParts.append('Part-{}'.format(i+1))
        p.BaseSolidExtrude(sketch=s, depth=layerD)
        s.unsetPrimaryObject()
        session.viewports['Viewport: 1'].setValues(displayedObject=p)

        del mdb.models['Model-1'].sketches['__profile__']
    except BaseException as e:
        print(e)
        pass


#plate
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(0,0), point2=block)
p = mdb.models['Model-1'].Part(name='Part-plate', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-plate']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
##
################
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
##########################################################
mc = []
for i,pn in enumerate(cutParts):
    p = mdb.models['Model-1'].parts[pn]
    a.Instance(name='{}-1'.format(pn), part=p, dependent=ON)
    mc.append(a.instances['{}-1'.format(pn)])
############
p = mdb.models['Model-1'].parts["Part-plate"]
a.Instance(name='Part-plate-1', part=p, dependent=ON)
# #
a = mdb.models['Model-1'].rootAssembly
a.InstanceFromBooleanCut(name='Part-base',
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Part-plate-1'],
    cuttingInstances=mc,
    originalInstances=DELETE)
# ####
cells= [a.instances["Part-base-1"]]
p = mdb.models['Model-1'].parts["Part-base"]
a.Instance(name='Part-base-2', part=p, dependent=ON)

###################1/4 to base model
a = mdb.models['Model-1'].rootAssembly
# a.translate(instanceList=('Part-base-2', ), vector=(block[0], 0.0, 0.0))
a.rotate(instanceList=('Part-base-2', ), axisPoint=(block[0], 0, 0.0),
    axisDirection=(0, 1, 0.0), angle=180.0)

#
# del mdb.models['Model-1'].parts["Part-base"]
a = mdb.models['Model-1'].rootAssembly
a.InstanceFromBooleanMerge(name='Part-a', instances=(
    a.instances['Part-base-1'], a.instances['Part-base-2'], ),
    originalInstances=DELETE, domain=GEOMETRY)
##########
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-a']
a.Instance(name='Part-a-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-a']
a.Instance(name='Part-a-2', part=p, dependent=ON)
#
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Part-a-2', ), axisPoint=(0.0, 0.0, 0.0),
    axisDirection=(1.0, 0.0, 0.0), angle=180.0)
a.translate(instanceList=('Part-a-2', ), vector=(0.0, 2.*block[0], 0.0))
#
a.InstanceFromBooleanMerge(name='Part-cell', instances=(
    a.instances['Part-a-1'], a.instances['Part-a-2'], ),
    originalInstances=DELETE, domain=GEOMETRY)
#######################
a = mdb.models['Model-1'].rootAssembly
a.LinearInstancePattern(instanceList=('Part-cell-1', ), direction1=(1.0, 0.0,
    0.0), direction2=(0.0, 1.0, 0.0), number1=int(pattern[0]/2.), number2=pattern[1], spacing1=block[0]*2.,
    spacing2=block[1]*2.)
###########
# ###########################################
a = mdb.models['Model-1'].rootAssembly
ins = a.instances.keys()
ins = [a.instances[k] for k in ins]
#####
a.InstanceFromBooleanMerge(name='Part-half', instances=ins, originalInstances=DELETE, domain=GEOMETRY)
#########
#############################################
p = mdb.models['Model-1'].parts['Part-half']
p.seedPart(size=meshSize, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-half']
f = p.faces
pickedRegions = f
p.setMeshControls(regions=pickedRegions, elemShape=TRI)
p = mdb.models['Model-1'].parts['Part-half']
p.generateMesh()


#################################################################
p = mdb.models['Model-1'].parts['Part-half']
p.PartFromMesh(name='Part-half1', copySets=True)
############
a = mdb.models['Model-1'].rootAssembly
# del a.features['Part-half-1']
#
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-half1']
a.Instance(name='Part-half-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-half1']
a.Instance(name='Part-half-2', part=p, dependent=ON)
########
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Part-half-2', ), axisPoint=(block[0]*2.*pattern[0]/2., 0.0, 0.0),
    axisDirection=(0., 1.0, 0.0), angle=180.0)
###########################
a = mdb.models['Model-1'].rootAssembly
a.InstanceFromBooleanMerge(name='Part-model', instances=(
    a.instances['Part-half-1'], a.instances['Part-half-2'], ), mergeNodes=ALL,
    nodeMergingTolerance=0.01, domain=MESH, originalInstances=DELETE)
############
############
#key********************************
p = mdb.models['Model-1'].parts['Part-model']
nodes = p.nodes
points = []
xs,ys ,zs = [],[],[]
for i, node in enumerate(nodes):
    coords = node.coordinates
    x, y,z = coords[0], coords[1],coords[2]
    points.append((x,y,z)) ######x,z,y
    xs.append(x)
    ys.append(y)
    zs.append(z)
#
xmin = min(xs)
xmax = max(xs)
total_length = abs(xmax - xmin)
R = total_length / (2 * math.pi)
circle_points = map_polyline_to_circle(points,R)
#
for i, node in enumerate(nodes):
    coords = node.coordinates
    ns = nodes[i:i+1]
    point = circle_points[i]
    p.editNode(nodes=ns, coordinate1=point[0],
        coordinate2=point[1], coordinate3=point[2])
###
#################################################################################
##################
D = 2.*R
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.Line(point1=(-D, H+0.1), point2=(D, H+0.1))
session.viewports['Viewport: 1'].view.setValues(width=228.74, height=74.4681,
    cameraPosition=(0.567771, -0.203095, 188.562), cameraTarget=(0.567771,
    -0.203095, 0))
p = mdb.models['Model-1'].Part(name='Part-top', dimensionality=THREE_D,
    type=DISCRETE_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['Part-top']
p.AnalyticRigidSurfExtrude(sketch=s, depth=D*2.0)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-top']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-top']
p.ReferencePoint(point=(0.0, H, 0))

#Part-bottom
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.Line(point1=(-D, -0.01), point2=(D, -0.01))
session.viewports['Viewport: 1'].view.setValues(width=228.74, height=74.4681,
    cameraPosition=(0.567771, -0.203095, 188.562), cameraTarget=(0.567771,
    -0.203095, 0))
p = mdb.models['Model-1'].Part(name='Part-bottom', dimensionality=THREE_D,
    type=DISCRETE_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['Part-bottom']
p.AnalyticRigidSurfExtrude(sketch=s, depth=D*2.)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-bottom']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bottom']
p.ReferencePoint(point=(0.0, 0.0, 0.0))
#########

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts["Part-top"]
a.Instance(name='Part-top-1', part=p, dependent=ON)
p = mdb.models['Model-1'].parts["Part-bottom"]
a.Instance(name='Part-bottom-1', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
# a.translate(instanceList=('Part-top-1', 'Part-bottom-1'), vector=(0.0, 0.0,
#     -D*2.))
###############
##
p = mdb.models['Model-1'].parts['Part-bottom']
p.seedPart(size=4.0, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()
p = mdb.models['Model-1'].parts['Part-top']
p.seedPart(size=4.0, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()
#material
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Density(table=((mat_data['density'], ), ))
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((mat_data['e'], mat_data['u']), ))
mdb.models['Model-1'].materials['Material-1'].Plastic(table=mat_data['plastic'])

# section
mdb.models['Model-1'].HomogeneousShellSection(name='Section-2',
    preIntegrate=OFF, material='Material-1', thicknessType=UNIFORM,
    thickness=Thickness_cell, thicknessField='', nodalThicknessField='',
    idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT,
    thicknessModulus=None, temperature=GRADIENT, useDensity=OFF,
    integrationRule=SIMPSON, numIntPts=5)
#
p = mdb.models['Model-1'].parts['Part-model']
e = p.elements
elements = e
region = regionToolset.Region(elements=elements)
p = mdb.models['Model-1'].parts['Part-model']
p.SectionAssignment(region=region, sectionName='Section-2', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)

#########

# step
a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['Part-top-1'].referencePoints
refPoints1=(r1[r1.keys()[-1]], )
a.Set(referencePoints=refPoints1, name='Set-top')

a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['Part-bottom-1'].referencePoints
refPoints1=(r1[r1.keys()[-1]], )
a.Set(referencePoints=refPoints1, name='Set-bottom')


mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial',
    timePeriod=step_time, massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0,
    7.8e-05, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))

############
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((
    0.2, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION,
    fraction=0.005, elasticSlipStiffness=None)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON,
    constraintEnforcementMethod=DEFAULT)
#
mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Step-1')
mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
    stepName='Step-1', useAllstar=ON)
mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
    stepName='Step-1', assignments=((GLOBAL, SELF, 'IntProp-1'), ))
#: The interaction "Int-1" has been created.

regionDef=mdb.models['Model-1'].rootAssembly.sets['Set-top']
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(variables=(
    'U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3', 'RF1', 'RF2', 'RF3', 'RM1', 'RM2',
    'RM3', 'ALLAE', 'ALLCD', 'ALLDMD', 'ALLFD', 'ALLIE', 'ALLKE', 'ALLPD',
    'ALLSE', 'ALLVD', 'ALLWK', 'ETOTAL'), frequency=1, region=regionDef,
    sectionPoints=DEFAULT, rebar=EXCLUDE)

#######merge node
p = mdb.models['Model-1'].parts['Part-model']
n = p.nodes
nodes = n
p.mergeNodes(nodes=nodes, tolerance=0.01)
##################
p = mdb.models['Model-1'].parts['Part-model']
nodes = p.nodes
top_nodes = None
bottom_nodes = None

ymin = 0
ymax = 2.*block[0]*pattern[1]

#
for i,node in enumerate(nodes):
    x,z = node.coordinates[0],node.coordinates[2]
    y = node.coordinates[1]
    if(abs(y-ymin)<=1e-5):
        if(bottom_nodes == None):
            bottom_nodes = nodes[i:i+1]
        else:
            bottom_nodes += nodes[i:i+1]
    if(abs(y-ymax)<=1e-5):
        if(top_nodes == None):
            top_nodes = nodes[i:i+1]
        else:
            top_nodes += nodes[i:i+1]

##
if(top_nodes !=None):
    p.Set(nodes=top_nodes, name='Set-node-top')

if(bottom_nodes !=None):
    p.Set(nodes=bottom_nodes, name='Set-node-bottom')
#########bc and load condition
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Part-model-1'].sets['Set-node-bottom']
mdb.models['Model-1'].EncastreBC(name='BC-bottom', createStepName='Initial',
    region=region, localCsys=None)

###
a = mdb.models['Model-1'].rootAssembly
region = a.sets['Set-bottom']
mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial',
    region=region, localCsys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
mdb.models['Model-1'].TabularAmplitude(name='Amp-1', timeSpan=STEP,
    smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (step_time, 1.0)))
a = mdb.models['Model-1'].rootAssembly
region = a.sets['Set-top']
mdb.models['Model-1'].DisplacementBC(name='BC-2', createStepName='Step-1',
    region=region, u1=0.0, u2=-displacement, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0,
    amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='',
    localCsys=None)
##########
###########################step/200
timeInterval = step_time/NNNN
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(
    timeInterval=timeInterval)
##############################################################
jobName = "Job-" + baseName
mdb.Job(name=jobName, model='Model-1')

#///
mdb.jobs[jobName].writeInput(consistencyChecking=OFF)