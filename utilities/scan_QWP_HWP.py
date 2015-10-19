#!/usr/bin/python
from subprocess import call, Popen, PIPE
import sys,math,os,time,re

def main():
    startHWP=0
    startQWP=0
    stepHWP=5
    stepQWP=5
    stopHWP=360
    stopQWP=360
    nSamples=10

    if sys.argv[1]=="default":
        print "You are going with default options!! Enjoy"
    elif len(sys.argv)==8:
        startHWP=sys.argv[1]
        stopHWP=sys.argv[2]
        stepHWP=sys.argv[3]
        startQWP=sys.argv[4]
        stopQWP=sys.argv[5]
        stepQWP=sys.argv[6]
        nSamples=sys.argv[7]
        print "Running with:"
        print "  HWP start point",startHWP
        print "  HWP start point",stopHWP
        print "  HWP start point",stepHWP
        print "  QWP start point",startQWP
        print "  QWP start point",stopQWP
        print "  QWP start point",stepQWP
        print "  # samples",nSamples
    else:
        print "Usage: ./scan_QWP_HWP.py [start postion for HWP <deg>] [stop position for HWP <deg>] [step size for HWP <deg>][start postion for QWP <deg>] [stop position for QWP <deg>] [step size for QWP <deg>] [# samples per step]"
        print " or: ./scan_QWP_HWP.py default"
        print "  default values are:"
        print "  HWP start point",startHWP
        print "  HWP start point",stopHWP
        print "  HWP start point",stepHWP
        print "  QWP start point",startQWP
        print "  QWP start point",stopQWP
        print "  QWP start point",stepQWP
        print "  # samples",nSamples
        sys.exit()    
        
    print "setting EPICS_CA_ADDR_LIST:"
    os.environ["EPICS_CA_ADDR_LIST"]="129.57.255.11 129.57.13.238 129.57.36.166 129.57.188.5 129.57.188.16 129.57.164.48 129.57.188.91"
    print os.environ["EPICS_CA_ADDR_LIST"]

    setPlates(startHWP,startQWP)

    #set step size
    moveHWP=abs(250*stepHWP)
    moveQWP=abs(250*stepQWP)
    call(["caput","COMPTON_PVAL_X_ao",moveHWP]) 
    call(["caput","COMPTON_PVAL_Z_ao",moveQWP]) 
    
    angleHWP=startHWP
    angleQWP=startQWP
    f=open("o_scan_HWP_QWP.txt","w")
    while (angleQWP<stopQWP):
        call(["caput","COMPTON_LMOVX_bo",moveQWP])
        angleQWP=angleQWP+stepQWP
        print "QWP now at ",angleQWP
        time.sleep(stepQWP*10)
        while (angleHWP<stopHWP):
            call(["caput","COMPTON_LMOVZ_bo",moveHWP])
            angleHWP=angleHWP+stepHWP
            print " HWP now at ",angleHWP
            time.sleep(stepHWP*10)
            S3a=[]
            S3b=[]
            for i in range(0,nSamples):
                S3a.append(caget("COMPTON_PSD_X",0))
                S3b.append(caget("COMPTON_PSD_Y",0))
                print "  !! reading",i," S3a S3b",S3a[-1],S3b[-1]

            (mA,dA)=detStat(S3a)
            (mB,dB)=detStat(S3b)
            if ( !math.isnan(dA) and !math.isnan(dB) ):
                print "  ~~ writing to file",angleQWP,angleHWP,mA,dA,mB,dB
                f.write(""+str(angleQWP)+" "+str(angleHWP)+" "+str(mA)+" "+str(dA)+" "+str(mB)+" "+str(dB)+"\n")
            else:
                print "  ?? problem with readings:",angleQWP,angleHWP,mA,dA,mB,dB
                
    f.close()
    ### could do an analysis module here!! FIXME
    
                
def caget(varName,varType):
    out=Popen(["caget","-t","-n",varName],stdout=PIPE)
    if varType==0:
        return float(out.stdout.read())
    elif varType==1:
        return int(out.stdout.read())
    elif varType==2:
        return out.stdout.read()
    else:
        print "caget: not sure what you want me to return: ", varName, varType
        sys.exit()

def setPlates(angleHWP, angleQWP):
    #setting to 0
    call(["caput","caput COMPTON_ORX_bo","1"])
    call(["caput","caput COMPTON_ORZ_bo","1"])
    stepHWP=caget("COMPTON_PVAL_X_ao",1)
    stepQWP=caget("COMPTON_PVAL_Z_ao",1)
    #250=1 deg
    moveHWP=abs(250*angleHWP)
    moveQWP=abs(250*angleQWP)
    call(["caput","COMPTON_PVAL_X_ao",moveHWP]) 
    call(["caput","COMPTON_PVAL_Z_ao",moveQWP]) 
    if angleHWP<0:
        call(["caput","COMPTON_RMOVX_bo",stepHWP])
    else:
        call(["caput","COMPTON_LMOVX_bo",stepHWP])

    if angleQWP<0:
        call(["caput","COMPTON_RMOVZ_bo",stepQWP])
    else:
        call(["caput","COMPTON_LMOVZ_bo",stepQWP])
    print "setting plates please wait ..."
    time.sleep(angleHWP*10)
    #set step size back to original
    call(["caput","COMPTON_PVAL_X_ao",stepHWP]) 
    call(["caput","COMPTON_PVAL_Z_ao",stepQWP]) 

    
        
def getStat(values):
    n = 0
    mean = 0.0
    M2 = 0.0
     
    for x in values:
        n = n + 1
        delta = x - mean
        mean = mean + delta/n
        M2 = M2 + delta*(x - mean)

    if n < 2:
        return (mean,float('nan'));
    else:
        return (mean,sqrt(M2 / (n - 1)))

    
if __name__ == '__main__':
    main()
                            




