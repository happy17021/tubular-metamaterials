#-*- coding:utf-8 -*-
import shutil
import sys
import os
import worker
#
sys.path.append('scripts')
#
import utils

#

# 
data = {}
data['temp'] = './temp'
data['result'] = os.path.join(data['temp'],'result')
data['data'] = os.path.join(data['temp'],'data')
data['images'] = './images'

# 
data['current'] = os.getcwd()

if(not os.path.exists(data['temp'])):
    os.mkdir(data['temp'])


if(not os.path.exists(data['result'])):
    os.mkdir(data['result'])

if (not os.path.exists(data['data'])):
    os.mkdir(data['data'])
# step 1: 
images = utils.getAllImages(data['images'])
# step 2: image 2 data
model_data = []
for i,image in enumerate(images):
    job_id,c = utils.image2data(image,data['data'])
    model_data.append((job_id,c))


# step 3:input
shutil.copy("scripts/genModel.py",os.path.join(data['temp'],"genModel.py"))

os.chdir(data['temp'])

all_jobs = []
for i,md in enumerate(model_data):
    try:
        Job_name = f"Job-{md[0]}.inp"
        if(os.path.exists(Job_name)):
            all_jobs.append(Job_name)
            continue
        print(f'th{i+1}model，name:{Job_name}')
        cmd = 'abaqus cae noGui=genModel.py -- {} {}'.format(md[0],md[1])
        os.system(cmd)
        all_jobs.append(Job_name)
    except BaseException as e:
        pass
