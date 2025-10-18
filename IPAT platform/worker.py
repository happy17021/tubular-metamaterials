# -*- coding:utf-8 -*-
import os
import shutil
import symtable
import sys
import threading
import time


class Worker(threading.Thread):
    def __init__(self,id = 0,job = "",args = {}):
        super(Worker, self).__init__()
        self.id = id
        self.job = job
        self.args = args


    def run(self):
        self.current_job = self.job[:-4]
        print(f"the{self.current_job}")
        cmd = ""
        cmd += f"call abaqus job=%s cpus={self.args['cpus']} ask=off" % self.current_job 
        os.system(cmd)
        time.sleep(60)
        while(True):
            if(os.path.exists(f"{self.current_job}.lck")):
                print(f'current job {self.current_job} is running!')
                time.sleep(10)
            else:
                print('Job finished !')
                self.getResult()
                break

        return

    def getResult(self):
        cmd = f'abaqus cae noGui=getOdb.py -- {self.current_job}'
        print(f'get odb: {cmd}')
        os.system(cmd)
        return