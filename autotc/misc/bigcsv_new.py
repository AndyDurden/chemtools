#!/bin/python
import shutil
import pprint
import subprocess
import csv
import sys
import os



class tocsv():
  def __init__(self, filename):
    self.filename=filename
    self.cols = []

  def casci(self):
    print("Getting casciE")
    self.get_casciE()
    print("Getting TDs")
    self.get_tds()
    print("Getting exchar")
    self.get_exchar()
    print("Writing file")
    self.writecsv()


  def writecsv(self):
    f = open(self.filename.split()[0]+".csv",'wb')
    wr = csv.writer(f)
    rows = [list(i) for i in zip(*self.cols)]
    for row in rows:
      wr.writerow(row)
    f.close()


  def get_tds(self):
    f = open(self.filename,'r')
    l = f.readline()
    rows = []
    while l:
      if l.split() != ["Singlet","state", "electronic", "transitions:"]: l = f.readline()
      else:
        print("found tds")
        f.readline() # the newline
        l = f.readline() # headers
        rows.append(l.split()[1:6])
        rows.append(["", "", "", "",""])
        f.readline();l=f.readline() # loads first row
        while l.split()[0] != "2":
          rows.append(l.split()[3:])
          l = f.readline()
    f.close()
    for col in [list(i) for i in zip(*rows)]:
      self.cols.append(col)

  def get_casciE(self):
    f = open(self.filename,'r')
    l = f.readline()
    rows = []
    while l:
      if l.split()!=['Root', 'Mult.', 'Total', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(eV)', 'Ex.', 'Energy', '(nm)']: l=f.readline()
      else:
        print("found Es")
        l = l.split() #headers
        rows.append([l[2]+l[3]+l[4], l[8]+l[9]+l[10]])
        f.readline();l=f.readline().split()
        rows.append([l[2], "0.0"])
        l=f.readline()
        while l!="\n":
          rows.append([l.split()[2],l.split()[3]])
          l = f.readline()
    f.close()
    for col in [list(i) for i in zip(*rows)]:
      self.cols.append(col)


  def get_exchar(self):
    f = open(self.filename,'r')
    l = f.readline().split()
    col = []
    while l:
      if l[0]+l[4]=="Transitiontrace":
	print("found exchar")
	col.append(l[4])
	f.readline();l=f.readline().split()
	while l[0]!="2":
	  col.append(l[3])
      else: l=f.readline().split()
    f.close()
    if len(col)>0:
      self.cols.append(col)
    else: print("Excitation Character not found.")

 

newcsv = True
jdir= "/home/adurden/jobs/dft_pno_test/"
def bigsheet(jdir):
  csvmat = [] # [col][row] of csv file 
  dirlist = next(os.walk(jdir))[1] # list of jobs generated by autojob
  # We are scanning parameters X and Y, and the directory name
  # for each job is Xval_Yval, so we make a dictionary (data) 
  # with keys [Xval], that stores dictionaries with keys [Yval], that stores the data.
  # So basically a 2D array, but with the directory strings as keys instead of integers.
  data = {}
  lencol = 0 
  for d in dirlist:
    data[d.split("_")[0]]={}
    for d1 in dirlist:
      data[d.split("_")[0]][d1.split("_")[1]]={}
  for subdir in dirlist:
    mat = []
    if ( (not os.path.isfile(jdir+subdir+"/ag4/ag4_"+subdir+".out")) or newcsv): 
      a = tocsv(jdir+subdir+"/ag4/ag4_"+subdir+".out")
      a.casci()
    f = open(jdir+subdir+"/ag4/ag4_"+subdir+".out.csv",'r')
    print(jdir+subdir+"/ag4/ag4_"+subdir+".out.csv")
    r = csv.reader(f, delimiter=",")
    for row in r: mat.append(row)

    print("mat:\n"+str(mat)+"\n")

    #print(map(list,zip(*mat)))
    csvmat = map(list,zip(*mat))
    if len(csvmat)>1:
      lencol = max(lencol, len(csvmat[0])) # failed jobs have column length 2 instead of nstates+1
    print("csvmat:\n"+str(csvmat)+"\n")
    f.close()
    #for i in range(1,lencol): data[subdir.split("_")[0]][subdir.split("_")[1]][i]={}
    for j in range(1,lencol):
     if len(csvmat)<6: pass # something screwy with csv if <6, skip.
     else: 
      data[subdir.split("_")[0]][subdir.split("_")[1]][j]=\
        {'e':csvmat[0][j], 'ex':csvmat[1][j], 'tx':csvmat[2][j],
        'ty':csvmat[3][j], 'tz':csvmat[4][j], 'tmag':csvmat[5][j],
        'tosc':csvmat[6][j], 'exc':None}
      if len(csvmat)>7:
        data[subdir.split("_")[0]][subdir.split("_")[1]][j]['exc']=csvmat[7][j]
  # I know this is very unelegant but it works:
  #lim=0
  #while len(bigdata[lim])<2: lim+=1
  #for i in range(1,len(bigdata[lim][0])):
  #  data[d.split("_")[0]][d1.split("_")[1]][i]={}
  #for i in range(0,len(dirlist)):
  #  for j in range(1,len(bigdata[lim][0])):
  #   if len(bigdata[i])<6: pass # something screwy with csv if <6, skip.
  #   else: 
  #    data[dirlist[i].split("_")[0]][dirlist[i].split("_")[1]][j]=\
  #      {'e':bigdata[i][0][j], 'ex':bigdata[i][1][j], 'tx':bigdata[i][2][j],
  #      'ty':bigdata[i][3][j], 'tz':bigdata[i][4][j], 'tmag':bigdata[i][5][j],
  #      'tosc':bigdata[i][6][j], 'exc':None}
  #    if len(bigdata[i])>7:
  #      data[dirlist[i].split("_")[0]][dirlist[i].split("_")[1]][j]['exc']=bigdata[i][7][j]
  #import pdb; pdb.set_trace() 
  print("Data structure done")
      

  # 'field' and 'iters' are just the two parameters scanned over
  #field = sorted(data.keys(), key=float)
  field = sorted(data.keys(), key=str.lower)
  iters = sorted(data[field[0]].keys(), key=float)

  csvdir = "/home/adurden/jobs/autotest/csv/dft_pno/"
  if not os.path.exists(csvdir): os.makedirs(csvdir)
  if not os.path.exists(csvdir+"e/"): os.makedirs(csvdir+"e/") # Eng (in au) of state N
  if not os.path.exists(csvdir+"ex/"): os.makedirs(csvdir+"ex/") # Ex. eng (in au) of state N
  if not os.path.exists(csvdir+"e_avg/"): os.makedirs(csvdir+"e_avg/") # Average of N states
  if not os.path.exists(csvdir+"t/"): os.makedirs(csvdir+"t/") # Transition dipole vector of state N
  if not os.path.exists(csvdir+"exc/"): os.makedirs(csvdir+"exc/") # Excitation character of state N
  for i in range(1,lencol):
    for jt in ['e','ex','exc','e_avg']:
      f = open(csvdir+jt+"/"+str(i)+".csv",'wb')
      wr = csv.writer(f)
      wr.writerow(['']+iters) #col headers
      for fi in field:
        row = [fi] #row headers
        for it in iters:
         try:
          if jt=='e_avg':
            avg = 0
            for i2 in range(1,i+1): avg+=float(data[fi][it][i2]['e'])
            avg = avg/i
            row.append(avg)
          else: row.append(data[fi][it][i][jt])
         except KeyError: row.append("")
        wr.writerow(row)
      f.close()

    f = open(csvdir+"t/"+str(i)+".csv",'wb')
    wr = csv.writer(f)
    row = ['']
    for it in iters: row=row+['Tx','Ty','Tz','|T|','Osc.','']
    wr.writerow(row)
    row = ['']
    for it in iters: row=row+[it,"","","","",""]
    wr.writerow(row)
    for fi in field:
      row = [fi] 
      for it in iters:
       try:
        row.append(data[fi][it][i]['tx'])
        row.append(data[fi][it][i]['ty'])
        row.append(data[fi][it][i]['ty'])
        row.append(data[fi][it][i]['tmag'])
        row.append(data[fi][it][i]['tosc'])
        row.append('')
       except KeyError: row=row+['','','','','','']
      wr.writerow(row)
    f.close()

#dirlist = next(os.walk(jdir))[1]
#for j in dirlist:
#  subprocess.call(["/home/adurden/Programs/bin/out2csv", jdir+j+"/ag4/ag4_"+j+".out"])

bigsheet(jdir)


