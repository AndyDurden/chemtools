#!/bin/python


import os
from maketest import *
import ConfigParser
import ast
print("\n")

settingsini = ConfigParser.ConfigParser()
settingsini.readfp(open("../settings.ini"))
settings = {}
for i in settingsini.items('config'):
  settings[i[0]] = i[1]

scheduler_mode = ast.literal_eval(settings["scheduler_mode"])
runjobs = ast.literal_eval(settings["runjobs"])
xyzdir = settings["xyzdir"]
jobdir = settings["jobdir"]
templatedir = settings["templatedir"]
if ("alt_template_name" not in settings.keys()) or (settings["alt_template_name"] == "None"):
  alt_template_name = None
else:
  alt_template_name = settings["alt_template_name"]



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
xyz_params.readfp(open("../xyz_params.ini"))
job_params = ConfigParser.ConfigParser()
job_params.readfp(open("../job_params.ini"))

for xyz in xyz_params.sections():
  for jobtype in job_params.sections():
    makefiles(xyz, jobtype, xyz_params, job_params, xyzdir, jobdir, templatedir)
    if runjobs:
      if scheduler_mode:
        print("Adding to queue: "+jobdir+"/"+xyz.split(".")[0]+"/"+jobtype+"/", xyz.split(".")[0]+"_"+jobtype+"\n")
        f = open('/proc/'+sched_pid+'/fd/0','a')
        if alt_template_name:
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
        f.close()
        #runjob(jobdir+"/"+jobtype+"/"+xyz.split(".")[0]+"/", xyz.split(".")[0]+"_"+jobtype)
        print("Added!\n\n")
      else: pass
        










