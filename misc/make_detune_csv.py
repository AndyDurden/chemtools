#!/bin/python

# Config variables here:
# ========================

# Use True on computers that aren't already using a scheduler, False if on a cluster (remember to use appropriate job templates)
scheduler_mode = True
# True to generate job input and run the jobs, False to only generate job input
runjobs = True

# Directory Settings:
# include leading slash on dirs

# All your .xyz files should go here:
xyzdir="/home/adurden/xstream/templates/xyz/"
# Terachem input templates go here. 
templatedir="/home/adurden/xstream/templates/in/"
# Directory where job input is written and executed
jobdir="/home/adurden/jobs/tdci_floquet_detune/"

# For jobs which use a different script to run than the default scripts set in the scheduler.
# alt_templates MUST run on the proper node (gpu) when the string "tempgpunum"
# is replaced by the node (gpu) number or bad things will happen. 
#alt_template_name = None
alt_template_name = "template_denali_rabi.job"

# End Config
# ========================

import os
from maketest import *
import ConfigParser
print("\n")

if scheduler_mode:
  # make sure scheduler is running
  if not os.path.exists("/tmp/schedaemon.pid"):
    print("scheduler daemon pid file not found\n")
    p = subprocess.Popen( "python scheduler.py start" , shell=True)
    p.wait()
    print("started the scheduler daemon")
  # get scheduler pid
  f = open('/tmp/schedaemon.pid','r')
  sched_pid = f.read()
  sched_pid = sched_pid.strip()
  f.close(); del f
  print("scheduler pid: "+sched_pid+"\n")



xyz_params = ConfigParser.ConfigParser()
xyz_params.readfp(open("xyz_params.ini"))
job_params = ConfigParser.ConfigParser()
job_params.readfp(open("job_params.ini"))

numA = {}
numT = {}
tdciA = {}
tdciT = {}



for xyz in xyz_params.sections():
  for jobtype in job_params.sections():
    f = open(jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/rabi_analysis.csv", 'r')
    f.readline();f.readline()
    a = f.readline()
    detune, step = jobtype.split("_")
    detune = float(detune.replace("n","-"))
    step = int(step)
    numA[ (detune, step) ] = float(a.split(",")[1])
    numT[ (detune, step) ] = float(a.split(",")[2])
    a = f.readline()
    tdciA[(detune,step)] = float(a.split(",")[1])
    tdciT[(detune,step)] = float(a.split(",")[2])
    f.close()
    del f
    """
    makefiles(xyz, jobtype, xyz_params, job_params, xyzdir, jobdir, templatedir)
    if runjobs == 1:
      if scheduler_mode:
        print("Adding to queue: "+jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/", xyz.split(".")[0]+"_"+jobtype+"\n")
        f = open('/proc/'+sched_pid+'/fd/0','a')
        if alt_template_name:
          #f.write( jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/,"+xyz.split(".")[0]+"_"+jobtype+","+alt_template_name+"\n")
          jobdict = {}
          jobdict["name"] = xyz.split(".")[0]+"_"+jobtype
          jobdict["path"] = jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/"
          jobdict["alt_template"] = alt_template_name
          f.write(str(jobdict)+"\n")
        else:
          jobdict = {}
          jobdict["name"] = xyz.split(".")[0]+"_"+jobtype
          jobdict["path"] = jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/"
          f.write(str(jobdict)+"\n")
          #f.write( jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/,"+xyz.split(".")[0]+"_"+jobtype+"\n")
        f.close()
        #runjob(jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/", xyz.split(".")[0]+"_"+jobtype)
        print("Added!\n\n")
      else: pass
    """ 

detunes = list(set([x[0] for x in numA.keys()]))
detunes.sort()
steps = list(set([x[1] for x in numA.keys()]))
steps.sort()

nA = open("numA.csv",'w')
nT = open("numT.csv",'w')
tA = open("tdciA.csv",'w')
tT = open("tdciT.csv",'w')
fs = [nA, nT, tA, tT]
for f in fs:
  f.write(",")
  for step in steps: f.write(str(step)+",")
  f.write("\n")
for detune in detunes:
  for f in fs:
    f.write(str(detune)+",")
  for step in steps:
    nA.write(str(numA[(detune,step)])+",")
    nT.write(str(numT[(detune,step)])+",")
    tA.write(str(tdciA[(detune,step)])+",")
    tT.write(str(tdciT[(detune,step)])+",")
  for f in fs:
    f.write("\n")
    

for f in fs:
  f.close()






