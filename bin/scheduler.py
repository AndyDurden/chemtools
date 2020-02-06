#!/bin/python

# Daemonizer from Sander Marechal
# http://web.archive.org/web/20131017130434/http://www.jejik.com/articles/2007/02/a_simple_unix_linux_daemon_in_python/

import sys, os, shutil, subprocess, time, atexit, ast
from daemon import Daemon
import signal

class schedaemon(Daemon):
  def __init__(self, pidfile):
    self.stdin = "/tmp/schedaemon_stdin"
    self.stdout = "/tmp/schedaemon_stdout"
    self.stderr = "/tmp/schedaemon_stderr"
    self.logdir = "/home/adurden/autotc/log/"
    self.logf = self.logdir+"log"
    self.pidfile = pidfile
    if not os.path.exists(self.stdin): open(self.stdin,'a').close()
    self.queue = []
    self.last_jobid = 0
    self.templatedir = "/home/adurden/autotc/templates/in/"
    self.nnodes = 3
    self.node_template = []
    self.node_proc = []
    self.node_job  = []
    for n in range(0,self.nnodes):
      (self.node_template).append("template_corecm.job")
      (self.node_proc).append(None)
      (self.node_job).append('')
    self.flushfiles()

  def flushfiles(self):
    f = open(self.stdin,'w');f.write("");f.close();del f
    f = open(self.stdout,'w');f.write("");f.close();del f
    f = open(self.stderr,'w');f.write("");f.close();del f
    i=1
    while( os.path.isfile(self.logf+str(i))): i+=1
    if( os.path.isfile(self.logf)): os.rename(self.logf, self.logf+str(i))
    f = open(self.logf,'w');f.write("");f.close();del f
    return 0

  def logprint(self,string):
    l = open(self.logf,'a')
    out = "["+time.asctime()+"] "+ string + "\n"
    l.write(out);l.close;del l

  def exec_job(self, qitem, node):
    if os.path.isfile(qitem["path"]+qitem["name"]+".job"): os.remove(qitem["path"]+qitem["name"]+".job")
    # 3rd string is alternate template
    if "alt_template" in qitem.keys():
      if os.path.isfile(self.templatedir+qitem["alt_template"]): shutil.copy( self.templatedir+qitem["alt_template"], qitem["path"]+qitem["name"]+".job") 
      else:
        l = open(self.logf,'a');l.write("["+time.asctime()+"] Couldn't find alternate template for job, skipping: "+str(qitem)+"\n" );l.close();del l
        return subprocess.Popen('exit 1', shell = True)
    else: 
      l = open(self.logf,'a');l.write("["+time.asctime()+"] Starting job with default template on node "+str(node)+": "+str(qitem)+"\n" );l.close();del l
      shutil.copy( self.templatedir+self.node_template[node], qitem["path"]+qitem["name"]+".job")
    p = subprocess.Popen( 'sed -i -e "s+temppath+'+qitem["path"]+'+g" '+qitem["path"]+qitem["name"]+'.job', shell=True)
    p.wait()
    p = subprocess.Popen( 'sed -i -e "s+tempname+'+qitem["name"]+'+g" '+qitem["path"]+qitem["name"]+'.job', shell=True)
    p.wait()
    # alt_templates MUST ensure their job runs on the appropriate node by replacing the text "tempgpunum" with the node number
    p = subprocess.Popen( 'sed -i -e "s+tempgpunum+'+str(node)+'+g" '+qitem["path"]+qitem["name"]+'.job', shell=True)
    p.wait()
    return subprocess.Popen( 'bash '+qitem["path"]+qitem["name"]+".job", shell=True,preexec_fn=os.setsid)


  # Put any triggers for restarting a job from output here.
  def check_restart_fail(self,qitem):
    f = open(qitem["path"]+qitem["name"]+".out",'r')
    line = f.readline()
    while line:
      # spin purification failure, not always reproducible in terachem.
      if line.split()[0:3] == ["Diagonalization","also","failed"]:
        return True
      line = f.readline()
    return False 


  # main daemon loop
  def run(self):
    while True:
      try:
        # Add jobs in stdin to queue
        line=sys.stdin.readline().rstrip()
        while line:
          self.logprint("Read from stdin: "+line)
          if line.split()[0] == "!command": # handle commands from stdin
            if len(line.split())<2: self.logprint("Empty command received, ignoring.")
            else: self.handle_command(line.split()[1:])

          # Add jobs directly added to stdin without command
          elif len(line)>0:
            if line[0] == "{":
              jobdict = ast.literal_eval(line.strip())
              if (("path" not in jobdict.keys()) or ("name" not in jobdict.keys())):
                l = open(self.logf,'a');l.write("["+time.asctime()+"] Badly formatted job dictionary, discarding: "+line+"\n" );l.close();del l
              elif (not (os.path.exists(jobdict["path"]))):
                l = open(self.logf,'a');l.write("["+time.asctime()+"] Job dictionary properly formatted, but path does not exist. Discarding: "+line+"\n" );l.close();del l
              else: # Well formatted job dictionary
                self.last_jobid += 1
                jobdict["jobid"] = self.last_jobid
                self.queue.append( jobdict )          
                self.logprint("Added to Queue: "+str(self.queue[-1]))
          else: l = open(self.logf,'a');l.write("["+time.asctime()+"] Improperly formatted line in stdin was ignored: "+str(line)+"\n" );l.close();del l
          line=sys.stdin.readline().rstrip()

        # Check if nodes are open, 
        for i in range(0,self.nnodes):
          if self.node_proc[i] == "Reserved": continue # skip reserved gpus
          if self.node_job[i] != '': # if empty, no job has run yet (since init) or last job was deleted
            if self.node_proc[i].poll() != None: # returns None if still running
              # Check if current job failed spin purification, if so add to end of queue to be redone
              if self.check_restart_fail( self.node_job[i] ):
                self.logprint("Failure warranting re-run, resubmitting: "+str(self.node_job[i]))
                self.queue.append( self.node_job[i])
              else:
                self.logprint("Job exited with code: "+str(self.node_proc[i].poll())+". Job: "+str(self.node_job[i]))
              if "reserve_after_completion" in self.node_job[i]:
                if self.node_job[i]["reserve_after_completion"]:
                  self.node_job[i] = ''
                  self.node_proc[i] = "Reserved"
                  self.logprint("Reserved node "+str(i)+" after job was completed.")
                  continue
              # Run next job
              if len(self.queue) > 0:
                self.node_job[i] = self.queue.pop(0)
                self.logprint("Starting next job: "+str(self.node_job[i])+" on node "+str(i)) 
                self.node_proc[i] = self.exec_job( self.node_job[i], i )
                self.logprint("Started job: "+str(self.node_job[i])+"on node "+str(i)+", PID: "+str(self.node_proc[i].pid))
              else:
                self.node_job[i] = '' # Queue is empty and job is finished, so empty node.
                self.node_proc[i] = None
                time.sleep(0.5) # nap while queue is empty
            else: time.sleep(0.1) # nap while job finishes
          else: # Run first job
            if len(self.queue) > 0:
              self.node_job[i] = self.queue.pop(0)
              self.logprint("Starting job on empty node: "+str(self.node_job[i])+" on node "+str(i))
              self.node_proc[i] = self.exec_job( self.node_job[i], i )
              self.logprint("Started job: "+str(self.node_job[i])+" on node "+str(i)+", PID: "+str(self.node_proc[i].pid))
            else: time.sleep(0.5) # nap while queue is empty
        time.sleep(0.1) # nap at end of main loop
      except Exception as err:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        l = open(self.logf,'a');l.write("["+time.asctime()+"] Big error: "+str(err)+"\n  Extra info: "+str((exc_type,exc_obj,exc_tb))+"\n");l.close();del l
        time.sleep(10) # sleep 10sec so hard drive doesnt fill up (as fast) if error loops


  def handle_command(self, command):
    if command[0] == "del":
      if len(command)<2: self.logprint("Delete request had no jobid, ignoring.")
      else: self.delete_job(command[1:])
    elif command[0] == "swap":
      if len(command)<3: self.logprint("Swap request didn't have two arguments, ignoring.")
      else: self.swap_jobs(int(command[1]),int(command[2]))
    elif command[0] == "jobinfo":
      if len(command)<2: self.logprint("Job Info request had no jobid, ignoring.")
      else: self.jobinfo(int(command[1]))
    elif command[0] == "reserve":
      if len(command)<2: self.logprint("Reservation request had no gpu number, ignoring.")
      else: self.reserve_node(command[1:])
    elif command[0] == "unreserve":
      if len(command)<2: self.logprint("Unreservation request had no gpu number, ignoring.")
      else: self.node_proc[int(command[1])] = None
      self.logprint("Node "+command[1]+" was unreserved by request.")
    elif command[0] == "logprintq":
      l = open(self.logf,'a');l.write("["+time.asctime()+"] Queue print requested:\n")
      for i in range(0,self.nnodes): l.write("R on node "+str(i)+": "+str(self.node_job[i])+"\n")
      for q in self.queue: l.write(str(q)+"\n")
      l.close();del l
    elif command[0] == "printq":
      self.logprint("Updating queuefile")
      f = open(self.logdir+"queue",'w')
      f.write("Queue as of "+time.asctime()+"\n")
      for i in range(0,self.nnodes):
        if self.node_proc[i] == "Reserved":
          f.write("R on node "+str(i)+": Reserved\n")
        else:
          f.write("R on node "+str(i)+": "+str(self.node_job[i])+"\n")
        f.flush;os.fsync(f.fileno())
      for q in self.queue:
        f.write(str(q)+"\n")
        f.flush;os.fsync(f.fileno())
      f.write("::ENDMESSAGE::")
      f.flush;os.fsync(f.fileno())
      f.close();del f

  def reserve_node(self, args):
    n = int(args[0])
    force = False
    if len(args) > 1:
      if args[1] == "force":
        force = True
    self.logprint("Reservation request for node "+str(n)+". force="+str(force))
    if force and self.node_proc[n] != None:
      self.logprint("Killing job on node "+str(n)+" to comply with reservation request."+str(self.node_job[n]))
      self.delete_job(self.node_job[n]['jobid'])
      self.logprint("Sending killed job to back of queue:"+str(self.node_job[n]))
      self.queue.append(self.node_job[n])
      self.node_job[n] = ''
      self.node_proc[n] = "Reserved"
      self.logprint("Set node "+str(n)+" to Reserved.")
    elif self.node_proc[n] != None: # How do we actually wait until the job is done?
      self.logprint("node "+str(n)+" is currently used by job :"+str(self.node_job[n])+". Waiting until job is finished to reserve node.")
      self.node_job[n]['reserve_after_completion'] = True
    else:
      self.node_job[n] = ''
      self.node_proc[n] = "Reserved"
      self.logprint("Set node "+str(n)+" to Reserved.")
    return None

  def swap_jobs(self, i,j):
    self.logprint("swapping jobs "+str(i)+" and "+str(j))
    self.logprint("Queue before swap:")
    l = open(self.logf,'a')
    for ii in range(0,self.nnodes): l.write("R on node "+str(ii)+": "+str(self.node_job[ii])+"\n")
    for q in self.queue: l.write(str(q)+"\n")
    l.close();del l
    il, jl = None, None
    for ii in range(0,len(self.queue)):
      if int(self.queue[ii]['jobid']) == i: il = ii
      elif int(self.queue[ii]['jobid']) == j: jl = ii
    if (il == None) or (ij == None):
      self.logprint("Could not find specified jobids. No swap performed.")
      return None
    else:
      tempitem = self.queue[il]
      self.queue[il] = self.queue[jl]
      self.queue[jl] = tempitem
      self.logprint("Queue after swap:")
      l = open(self.logf,'a')
      for ii in range(0,self.nnodes): l.write("R on node "+str(ii)+": "+str(self.node_job[ii])+"\n")
      for q in self.queue: l.write(str(q)+"\n")
      l.close();del l
      return None

  def jobinfo(self, jobid):
    return None

  def delete_job(self, jobid):
    if type(jobid) == type([]): jid = int(jobid[0])
    else: jid = int(jobid)
    for q in self.queue: # Remove jobid from queue
      if q["jobid"] == jid:
        self.logprint("Removing queue item: "+str(q))
        self.queue.remove(q)
    for n in range(0,self.nnodes): # If jobid is running on a node, kill it.
      if type(self.node_job[n]) == type({}): # if node is unused, 'job' is empty string
        if self.node_job[n]["jobid"] == jid:
          self.logprint("Sending SIGTERM to process on node "+str(n)+": "+str(self.node_job[n]))
          #(self.node_proc[n]).terminate() # send SIGTERM (request kill)
          os.killpg(os.getpgid( (self.node_proc[n]).pid), signal.SIGTERM)
          time.sleep(5.0) # wait 5 secs for process to clean up
          if (self.node_proc[n]).poll() == None: # didn't die
            self.logprint("Process did not respond to SIGTERM within 5 seconds. Sending SIGKILL to process on node "+str(n)+": "+str(self.node_job[n]))
            #(self.node_proc[n]).kill() # send SIGKILL (kernel force close/kill -9)
            os.killpg(os.getpgid( (self.node_proc[n]).pid), signal.SIGKILL)
            time.sleep(5.0) # wait 5 secs for os to work its magic? idk if this is necessary
            
            if (self.node_proc[n]).poll() == None: # didn't die from SIGKILL
              self.logprint("SIGKILL didn't work?! Can't get rid of process with pid "+str((self.node_proc[n]).pid)+". Sorry but you'll need to look into this yourself.")
            else: # process ended from SIGKILL
              self.logprint("SIGKILL successful, process killed.")
          else: # process has ended from SIGTERM
            self.logprint("SIGTERM successful, process killed.")
          # Don't need to overwrite node_proc[n] or node_job[n], daemon should
          # detect that the process has ended on the next loop.

 
if __name__ == "__main__":
        daemon = schedaemon('/tmp/schedaemon.pid')
        if len(sys.argv) == 2:
                if 'start' == sys.argv[1]:
                        daemon.start()
                elif 'stop' == sys.argv[1]:
                        daemon.stop()
                elif 'restart' == sys.argv[1]:
                        daemon.restart()
                else:
                        print "Unknown command"
                        sys.exit(2)
                sys.exit(0)
        else:
                print "usage: %s start|stop|restart" % sys.argv[0]
                sys.exit(2)
