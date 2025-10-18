#-*- coding:utf-8 -*-
import glob
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
data['cpus'] = 2 #cpus
data['gpus'] = 0 #gpu
data['threads'] = 20 #CPU/2


#
data['current'] = os.getcwd()
data['result'] = os.path.join(data['temp'],'result')
if(not os.path.exists(data['temp'])):
    os.mkdir(data['temp'])

if(not os.path.exists(data['result'])):
    os.mkdir(data['result'])

#
all_jobs = glob.glob(os.path.join(data['temp'],"*.inp"))
all_jobs = [os.path.basename(i) for i in all_jobs]
# step 3:input
shutil.copy("scripts/getOdb.py",os.path.join(data['temp'],"getOdb.py"))
os.chdir(data['temp'])


cmd = "abaqus cae noGUI=getOdb.py"
os.system(cmd)

# sys.exit(0)
#
#
# 
threads = []
count = 1
#//////////////////////////////////
for i in range(data['threads']):
    job_name = all_jobs.pop(0)
    if(os.path.exists(job_name+'.odb')):
        continue
    thread = worker.Worker(count, job_name, data)
    thread.start()
    threads.append(thread)
    count += 1
    if(len(all_jobs) == 0):
        break

# 
while(True):
    for i,thread in enumerate(threads):
        if(thread.is_alive()):
            continue
        else:
            # 
            thread = threads.pop(i)
            #/////////////////////
            if(len(all_jobs) == 0):
                break
            job_name = all_jobs.pop(0)
            thread = worker.Worker(count, job_name, data)
            threads.append(thread)
            thread.start()
            count += 1
    if(len(threads) == 0):
        break
print("All Jobs have completed")
#

