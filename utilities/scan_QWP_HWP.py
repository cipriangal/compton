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

    if len(sys.argv)==2 and sys.argv[1]=="default":
        print "You are going with default options!! Enjoy"
    elif len(sys.argv)==8:
        startHWP=sys.argv[1]
        stopHWP=sys.argv[2]
        stepHWP=sys.argv[3]
        startQWP=sys.argv[4]
        stopQWP=sys.argv[5]
        stepQWP=sys.argv[6]
        nSamples=sys.argv[7]
    else:
        print "Usage: ./scan_QWP_HWP.py [start postion for HWP <deg>] [stop position for HWP <deg>] [step size for HWP <deg>][start postion for QWP <deg>] [stop position for QWP <deg>] [step size for QWP <deg>] [# samples per step]"
        print " or: ./scan_QWP_HWP.py default"
        print "  default values are:"
        print "  HWP start point",startHWP
        print "  HWP stop point",stopHWP
        print "  HWP step size",stepHWP
        print "  QWP start point",startQWP
        print "  QWP stop point",stopQWP
        print "  QWP step size",stepQWP
        print "  # samples",nSamples
        sys.exit()

    print "Running with:"
    print "  HWP start point",startHWP
    print "  HWP stop point",stopHWP
    print "  HWP step size",stepHWP
    print "  QWP start point",startQWP
    print "  QWP stop point",stopQWP
    print "  QWP step",stepQWP
    print "  # samples",nSamples

    print "setting EPICS_CA_ADDR_LIST:"
    os.environ["EPICS_CA_ADDR_LIST"]="129.57.255.11 129.57.13.238 129.57.36.166 129.57.188.5 129.57.188.16 129.57.164.48 129.57.188.91"
    print os.environ["EPICS_CA_ADDR_LIST"]

    angleQWP=startQWP
    f=open("o_scan_HWP_QWP.txt","w")
    while (angleQWP<=stopQWP):
        setPlate(startQWP,"QWP")
        angleQWP=angleQWP+stepQWP

        angleHWP=startHWP
        while (angleHWP<=stopHWP):
            setPlate(angleHWP,"HWP")
            angleHWP=angleHWP+stepHWP

            S3a=[]
            S3b=[]
            for i in range(0,nSamples):
                S3a.append(caget("COMPTON_PSD_X",0))
                S3b.append(caget("COMPTON_PSD_Y",0))
                print "  !! reading",i," S3a S3b",S3a[-1],S3b[-1]
                time.sleep(3)

            (mA,dA)=getStat(S3a)
            (mB,dB)=getStat(S3b)
            if ( not math.isnan(dA) and not math.isnan(dB) ):
                print "  ~~ writing to file",angleQWP,angleHWP,mA,dA,mB,dB
                f.write(""+str(angleQWP)+" "+str(angleHWP)+" "+str(mA)+" "+str(dA)+" "+str(mB)+" "+str(dB)+"\n")
            else:
                print "  ?? problem with readings:",angleQWP,angleHWP,mA,dA,mB,dB

    f.close()

    ## Set readback to default 5
    call(["caput","COMPTON_SURUGA_RBRATE_QW1","5"])
    call(["caput","COMPTON_SURUGA_RBRATE_QW2","5"])
    call(["caput","COMPTON_SURUGA_RBRATE_HW1","5"])

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

def setPlate(angle,plate):
    val=""
    if plate=="QWP":
        val="Z"
    else:
        val="X"

    setReadBack(plate)

    currentValue=caget("COMPTON_PVAL_"+val+"_ao",1)##FIXME which variables are we getting from???
    move = currentValue - 250*angle

    defaultStep=caget("COMPTON_PVAL_"+val+"_ao",1)

    call(["caput","COMPTON_PVAL_"+val+"_ao",str(move)])
    if move<0:
        call(["caput","COMPTON_RMOV"+val+"_bo","1"])
    else:
        call(["caput","COMPTON_LMOV"+val+"_bo","1"])
    print "setting plates please wait ...",angle*0.2," sec"
    time.sleep(angle*0.2)

    confirmedMove = False
    nLoop = 0
    while not confirmedMove:
        currentValue = caget("COMPTON_PVAL_"+val+"_ao",1)
        if currentValue == angle*250:
            confirmedMove = True
        nLoop++
        time.sleep(1)
        if nLoop>10:
            setPlate(angle,plate)
    print "Moved ",plate," to ", angle
    call(["caput","COMPTON_PVAL_"+val+"_ao",str(defaultStep)])

def setReadBack(plate):
    if plate == "QWP"
        call(["caput","COMPTON_SURUGA_RBRATE_QW1","1"])
        call(["caput","COMPTON_SURUGA_RBRATE_QW2","0"])
        call(["caput","COMPTON_SURUGA_RBRATE_HW1","0"])
    elif plate == "HWP":
        call(["caput","COMPTON_SURUGA_RBRATE_QW1","0"])
        call(["caput","COMPTON_SURUGA_RBRATE_QW2","0"])
        call(["caput","COMPTON_SURUGA_RBRATE_HW1","1"])
    else:
        print "I don't know what you want me to do with ", plate

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
        return (mean,math.sqrt(M2 / (n - 1)))


if __name__ == '__main__':
    main()





