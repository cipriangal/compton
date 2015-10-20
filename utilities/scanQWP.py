from subprocess import call, Popen, PIPE
import os
import time
import re

print "setting EPICS_CA_ADDR_LIST:"
os.environ["EPICS_CA_ADDR_LIST"]="129.57.255.11 129.57.13.238 129.57.36.166 129.57.188.5 129.57.188.16 129.57.164.48 129.57.188.91"
print os.environ["EPICS_CA_ADDR_LIST"]

#set step size
call(["caput","COMPTON_PVAL_X_ao","750"]) 
out=Popen(["caget","-t","-n","COMPTON_PVAL_X_ao"],stdout=PIPE)
stepSize=int(out.stdout.read())
startPos=42
print "Starting scan at ",startPos," deg and stepping with ",stepSize/250," deg"

nstep=360/(stepSize/250)+1
f=open("o_scan.txt","w")
for i in range(1,nstep):
    deg=(startPos-i*stepSize/250+360)%360
    print "step:",i," rotating to ",deg," deg"
    call(["caput","COMPTON_RMOVX_bo","1"])
    time.sleep(10)
    out=Popen(["caget","-t","-n","COMPTON_PW1R_S2_ca"],stdout=PIPE)
    pow2=float(out.stdout.read())
    out=Popen(["caget","-t","-n","COMPTON_PW1R_S1_ca"],stdout=PIPE)
    pow1=float(out.stdout.read())
    out=Popen(["caget","-t","-n","COMPTON_PW2R_S2_ca.VAL"],stdout=PIPE)
    norm=float(out.stdout.read())
    out=Popen(["caget","-t","-n","COMPTON_SERVO_PDRai"],stdout=PIPE)
    pdr=float(out.stdout.read())
    f.write(""+str(deg)+" "+str(pow1)+" "+str(pow2)+" "+str(norm)+" "+str(pdr)+"\n")
    print deg,pow1,pow2,norm,pdr
f.close()

call(["root","-l","-b","-q","analyzeScan.C"])




