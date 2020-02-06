#!/bin/python

import sys
import os
import glob
import subprocess
import time
import ast

schedpath = os.getenv('SCHEDAEMONPATH',"/home/adurden/autotc/bin/")

def get_sched_pid():
  if not os.path.exists("/tmp/schedaemon.pid"): return False
  else:
    f = open('/tmp/schedaemon.pid','r')
    pid = (f.read()).strip()
    f.close()
    return pid

def start_sched():
  pid = get_sched_pid()
  if pid:
    print("Schedaemon already started\n")
    return pid
  else:
    p = subprocess.Popen(["python",schedpath+"scheduler.py","start"])
    p.wait()
    return get_sched_pid()
    

def stop_sched():
  pid = get_sched_pid()
  if not pid:
    print("Schedaemon wasn't running (no pidfile)\n")
    return pid
  else:
    p = subprocess.Popen(["python",schedpath+"scheduler.py","stop"])
    p.wait()
    return get_sched_pid()


def restart_sched():
  pid = get_sched_pid()
  if pid:
    stop_sched()
    return start_sched()
  else:
    print("Schedaemon wasn't running (no pidfile). Starting\n")
    return start_sched()


def addjob(schedpid,argstring):
  jobdict = {}
  if argstring == "": # Create job from current directory
    jobdict["path"] = (os.getcwd())+"/"
    infiles = glob.glob('*.in')
    if len(infiles)==0:
      print("No input file found in current directory\n")
      return False
    if len(infiles)>1:
      print("Multiple input files found. Using: "+str(infiles[0])+"\n")
    jobdict["name"] = (infiles[0]).split(".")[0] 
  # Really difficult to format the dictionary properly on command line, so forget this.
  #elif argstring[0] == "{": # Assume its a proper literal
  #  jobdict = ast.literal_eval(argstring)  
  #  if "path" not in jobdict.keys:
  #    print("Input string didnt have a job path?\n")
  #    return False
  #  if "name" not in jobdict.keys:
  #    print("Input string didnt have a job name?\n")
  #    return False
  else: # assume argstring is alt_template
    jobdict["path"] = (os.getcwd())+"/"
    infiles = glob.glob('*.in')
    if len(infiles)==0:
      print("No input file found in current directory\n")
      return False
    if len(infiles)>1:
      print("Multiple input files found. Using: "+str(infiles[0])+"\n")
    jobdict["name"] = (infiles[0]).split(".")[0] 
    jobdict["alt_template"] = argstring
    

  f = open('/proc/'+schedpid+'/fd/0','a')
  f.write(str(jobdict)+"\n")
  f.close()
  return True
