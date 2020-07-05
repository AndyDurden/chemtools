import os
import shutil
import pprint
import subprocess

# This file has some miscillaneous functions, most importantly for
# setting up terachem job inputs
# there's also some deprecated functions which basically diff two terachem outputs


# Used this function in gputest where I wanted a static active space
# for different molecules, but now its unused
# count electrons in molecule from xyz file
# assuming all electrons are treated explicitly
def counte(xyz):
  ## All electrons explicit
  #d = {"H" : 1,  "LI":3, "B":5, "C":6, "N":7, "O":8, "F":9,
  #     "NA":11, "MG":12, "P":15, "S":16, "CL":17}

  # Lanl2dz ecp
  d = {"H" : 1,  "LI":3, "B":5, "C":6, "N":7, "O":8, "F":9,
        "NA":1, "MG":2, "AL":3, "SI":4, "P":5, "S":6, "CL":7,
        "K":9, "CA":10, "SC":11, "TI":12, "V":13, "CR":14, "MN":15,
        "FE":16, "CO":17, "NI":18, "CU":19, "ZN":12, "GA":3, "GE":4,
        "AS":4, "SE":6, "BR":7, "RB":9, "SR":10, "Y":11, "ZR":12, "NB":13,
        "MO":14,"TC":15,"RU":16,"RH":17,"PD":18,"AG":19, "CD":12, "IN":3,
        "SN":4, "SB":5, "TE":6, "I":7, "OS":16, "IR":17, "PT":18, "AU":19,
        "HG":12 }


  f = open(xyz,'r')
  f.readline();f.readline()
  e = 0
  cur = f.readline()
  while cur != "":
    a = cur.split()[0].upper()
    e += d[a]
    cur = f.readline()
  print(e)
  return e


# note: jobdir is the directory where all job subdirectories are stored
# jobpath is the subdirectory for this specific job
def makefiles(xyz, job, xyz_params, job_params, xyzdir, jobdir, templatedir):
  jobpath = jobdir+"/"+xyz.split(".")[0]+"/"+job+"/"
  jobname = xyz.split(".")[0]+"_"+job
  if not os.path.exists(jobpath): os.makedirs( jobpath )
  #shutil.copy( 'input_templateemplate.in' , jobpath+jobname+'.in')
  shutil.copy( xyzdir+xyz , jobpath+xyz)

  #input_template = open(templatedir+job+".in", 'r')
  input_template = open(templatedir+job_params.get(job,"template"), 'r')
  input_file = open(jobpath+jobname+".in", 'w')
  for line in input_template:
    if len(line.split())>0:
      if line.split()[0] == "coordinates":
        input_file.write("coordinates    "+xyz+"\n")
      else:
        matchany = False
        # substitute job specific parameters
        for option in job_params.options(job):
          if line.split()[0] == option:
            matchany = True
            # option value in job_params is option key in xyz_params
            input_file.write(option+"    "+xyz_params.get(xyz, (job_params.get(job, option)) )+"\n")
        if not matchany: input_file.write(line) # preserve line if not altered
    else: input_file.write(line) # preserve empty lines
  input_file.close()
  input_template.close()

  # make .job file
  shutil.copy( templatedir+"template_corecm_floqci.job", jobpath+jobname+".job")
  p = subprocess.Popen( 'sed -i -e "s+temppath+'+jobpath+'+g" '+jobpath+jobname+'.job', shell=True)
  p.wait()
  p = subprocess.Popen( 'sed -i -e "s+tempname+'+jobname+'+g" '+jobpath+jobname+'.job', shell=True)
  p.wait()

  return 0

def runjob(jobpath, jobname):
  #p = subprocess.Popen(  './tc.sh %s %s %s %s' % (str(jobdir),str(gpu),str(infile),str(outfile),)  , shell=True)
  p = subprocess.Popen( 'bash '+jobpath+jobname+".job", shell=True)
  p.wait()
  return 0


# deprecated stuff below here
# diff command that skips to final data in output file for hardware difference/reproducibility tests
def comparecis(out1, out2, logpath, jobtype):
  f1_ = open(out1, 'r')
  f2_ = open(out2, 'r')
  log = open(logpath, 'w')
  f1, f2 = [], []
  # This will break if an outfile happens to be larger than ram, but that never happens
  for line in f1_: f1.append(line) 
  for line in f2_: f2.append(line)
  f1ran = True
  f2ran = True
  # Make sure the outfiles are longer than 10 lines
  # so that checking if it finished properly doesnt
  # give index errors
  print("lens: "+ str(len(f1)) +"  "+str(len(f2)))
  if len(f1)<11:
    f1ran = False
    log.write(out1+"didnt run\n")
  if len(f2)<11:
    f2ran = False
    log.write(out2+"didnt run\n")
  if ( (not f1ran) or (not f2ran) ):
    log.close()
    return 1

  # Check if the job finished properly
  f1_finished = False
  f2_finished = False
  for l in f1[-10:]:
    if len(l.split())<2:pass
    elif l.split()[0]+l.split()[1] == "Jobfinished:": f1_finished = True
  for l in f2[-10:]:
    if len(l.split())<2:pass
    elif l.split()[0]+l.split()[1] == "Jobfinished:": f2_finished = True
  
  if ( (not f1_finished) or (not f2_finished)):
    if (not f1_finished): log.write(out1+" didn't finish\n")
    if (not f2_finished): log.write(out2+" didn't finish\n")
    log.close()
    return 1

  log.write("Both files finished properly\n")
  
  i1,j1 = 0,0
  i2,j2 = 0,0

  # This clause is all that needs to be changed to add more jobtypes
  # See getoutputbounds comments
  if jobtype == "cis":
    i1,j1 = getoutputbounds(f1,"RootTotalEnergy", "Totalprocessingtime:")
    i2,j2 = getoutputbounds(f2,"RootTotalEnergy", "Totalprocessingtime:")
  elif jobtype == "cisno":
    i1,j1 = getoutputbounds(f1,"Singletstate1", "TotalCASCITime")
    i2,j2 = getoutputbounds(f2,"Singletstate1", "TotalCASCITime")
  else:
    print("???")
    log.write("unknown job type???")
    log.close()
    return 1
  
  if f1[i1:j1] == f2[i2:j2]:
    log.write("All "+jobtype+" Outputs match perfectly! :D\n")
  else:
    log.write("Final "+jobtype+" Outputs have mismatch, printing lines that differ with "+out1+" first and "+out2+" second:\n\n")
    if j1-i1 != j2-i2:
      log.write("Something weird happened: Data output blocks of two output files have different numbers of lines??? j1-i1="+str(j1-i1)+" , j2-i2="+str(j2-i2)+" :\n\n")
    i,j = i1,i2
    while i < min(j1,j2):
      if f1[i] != f2[j]: log.write("line "+str(i)+": "+f1[i]+"line "+str(j)+": "+f2[j]+"\n")
      i+=1;j+=1
  print("Finished comparing: "+out1+" "+out2)
  log.close()
  return 0


# some good keys have 3 words and only appear once in an output file
# so startkey and endkey must be the flist.split()[0]+flist.split()[1]+flist.split()[2]
def getoutputbounds(flist, startkey, endkey):
  # Find start of output data
  i=0
  done = 0
  while done == 0:
    if len(flist[i].split())<3:i+=1
    elif flist[i].split()[0]+flist[i].split()[1]+flist[i].split()[2] == startkey: done=1
    else: i+=1
    if i>len(flist): done=1
    
  # Find end of output data
  j=i
  done = 0
  while done == 0:
    if len(flist[j].split())<3:j+=1
    elif flist[j].split()[0]+flist[j].split()[1]+flist[j].split()[2] == endkey: done=1
    else: j+=1
    if j>len(flist): done=1
  return i,j





