# GNU GENERAL PUBLIC LICENSE

#   **  NS_pocket  Copyright (C) 2022  Luca Gagliardi
#       Affiliation: Istituto Italiano di Tecnologia
#
# DESCRIPTION:
#   Interface for the NanoShaper pocket detection software with volume ranking of pockets
#   NanoShaper info: http://www.plosone.org/article/metrics/info%3Adoi%2F10.1371%2Fjournal.pone.0059744



import subprocess
import os
import numpy as np
import re
from datetime import datetime
import sys


runFolder = os.path.abspath('temp')+'/'
resultFolder = 'results/'

infoFileName = runFolder+'cavAtomsSerials.txt'
volumeFileName = runFolder+'cavitiesSize.txt'



def readPocketAtoms(atomsLines,resMap):
    '''
    From the indexes of the xyzr file, returns the desired portion of PQR. 
    Need as input the indexes and the whole protein residues info
    '''
    atomsLines = atomsLines.split()
    # print(atomsLines)
    subSet = [resMap[int(i)-1] for i in atomsLines]

    return subSet


def get_volume(lines,pnumber):
    # lines = volumeFile.readlines()
    # print(lines)
    volume = 0
    s =0
    for line in lines[1:]:
        # print(line)
        l = line.split()
        isPocket = not int(l[-1])
        # print(isPocket)
        if isPocket:
            # print(line)
            volume = float(l[-2])
            # print(lines[s + 7])
            # print(volume)
            if s==pnumber:
                # print('searched pocket OK')
                break
            s+=1
    if (volume==0):
        raise Exception("Cannot find volume")
    return volume

def get_proteinAtoms(filename):
    '''
    Use PQR to filter out undesired atoms such as H or in any case keep atoms meaningful to the SES
    '''
    try:
        # print(structure+'.pdb')
        inFile = open(filename,'r')
    except Exception:
        raise NameError("Cannot load PQR file")
    # try:
    #     # print(structure+'.pdb')
    #     _check = open(structure+'.pdb','r')
    # except Exception:
    #     raise NameError("Cannot load PDB file")
    comment =['#', 'CRYST[0-9]?']
    remark = ['REMARK']
    termination = ['TER', 'END', '\n']
    skip = comment+remark+termination
    skip = '(?:% s)' % '|'.join(skip)
    for line in inFile: 
        if(re.match(skip,line)): 
            pass 
        else:
            linegNOChain=re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            linegChain = re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            break

    if(linegChain):
        # print("PQR contains CHAIN_ID")
        isChainID=1                                                        #resID
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
    elif(linegNOChain):
        # print("PQR does NOT contain CHAIN_ID")
        isChainID =0                                        # resID
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
    else:
        raise NameError("Incorrect pqr file formatting")
    if(isChainID):
        resInd = 5
        chargeInd = 9
        rInd = 10
    else:
        resInd = 4
        chargeInd = 8
        rInd = 9
    nameInd = 3
    atomInd = 2
    coordInd = resInd +1
    
    inFile.seek(0)

    resMap=[]
    for line in inFile:
        if(re.match(skip,line)): 
            pass 
        else: 
            mline=re.match(matchPattern,line).groups()
            content = {'resName':mline[nameInd],'resNum':mline[resInd],'resAtom': mline[atomInd],'atomNumber':int(mline[1]),
                    'charge':float(mline[chargeInd]),'coord':list(map(float, mline[coordInd:coordInd+3])),'radius':float(mline[rInd])}
            resMap.append(content)

    
    # coords = np.array([p['coord'] for p in resMap])
    return resMap


def savePQR(atomsLines,resMap,ind,outFolder):
    '''
    Creates a PQR of the structure atoms contacted by the pocket
    '''
    filename = outFolder+'p'+str(ind+1)+"_atm.pqr"
    outFile = open(filename,'w')
    for ind in atomsLines:
        try:
            rChain = resMap[ind]['resChain']
        except KeyError:
            rChain = 'A'
        outFile.write("{:<6s}{:>5d} {:<5s}{:>3s} {:1s}{:>5s}   {:>8.3f} {:>8.3f} {:>8.3f} {:>8.4f} {:>8.4f}\n".format('ATOM',resMap[ind]['atomNumber'],
        resMap[ind]['resAtom'],resMap[ind]['resName'],rChain,resMap[ind]['resNum'],
        resMap[ind]['coord'][0],resMap[ind]['coord'][1],resMap[ind]['coord'][2],resMap[ind]['charge'],resMap[ind]['radius']))
    outFile.close()


# ===========================================================
# ===========================


def main():
    import getopt
    argv = sys.argv[1:]
    makeOFF=False
    try:
        opts, args = getopt.getopt(argv,"h",["off","help","boolMAP"])
    except getopt.GetoptError:
        print ('<<ERROR>> Uncorrect formatting of options or unavaible option')
        sys.exit(2)
    # print('options:',opts)
    # print('arguments:',args)
    for opt, arg in opts:
        if opt in ["-h","--help"]:
            print('Usage:\npython3 NS_pocketVol <structure_name> (in PQR format)')
            print("OPTIONS:")
            print("--off: stores pockets in off format")
            print("--boolMAP: a single _labelled.pqr file with dummy integers in the Charge field of the original PQR:\n\t0=not associated to any binding site.\n\t1-->10 belonging to pocket with rank given by the integer")
            input('\n')
            sys.exit()
        elif opt == '--off':
            makeOFF=True
        elif opt == '--boolMAP':
            makeMAP=True

    filename = args[0] #sys.argv[1:][0]
    match =re.match("(.+)\.pqr",filename)
    if(match):
        structureName = match.groups()[0]
        print(structureName)
        # input()
        pass
    else:
        print("ABORT: must provide a pqr file")
        exit()
    
    isFolder = os.path.isdir(resultFolder)
    if not isFolder:
        subprocess.run(['mkdir',resultFolder])

    start = datetime.now()
    # LOADING PROTEIN
    # print(filename)
    resMap = get_proteinAtoms(filename)
    # print(len(resMap))
    # print(proteinAtoms)
    protein_atoms = np.empty((0,4))
    for i in resMap:
        c = np.append(np.asarray(i['coord']),i['radius'])
        protein_atoms = np.vstack([protein_atoms,c])
    # NS XYZR FILE
    np.savetxt(runFolder+"NanoShaper_input.xyzr",protein_atoms,delimiter="\t",fmt='%.4f')

    # NS RUN

    try:
        out = subprocess.check_output(['./NanoShaper','template.prm'],cwd=runFolder)
    except subprocess.CalledProcessError as grepexc:                                                                                                   
        print ("NS error code", grepexc.returncode, grepexc.output)
        exit()

    #Create output folder 

    isFolder = os.path.isdir(resultFolder+structureName)
    
    if not isFolder:
        subprocess.run(['mkdir',resultFolder+structureName])
    
    outFolder = resultFolder+structureName+'/'
    print(outFolder)

    infoFile = open(infoFileName,'r')
    infoLines = infoFile.readlines()
    volumeFile = open(volumeFileName,'r')
    volumeFileLines = volumeFile.readlines()

    # RANKING 
    vol=[]
    for r in range(len(infoLines)):
        vol.append(get_volume(volumeFileLines,pnumber = r))
    # print(vol)
    sortedInd=np.argsort(vol)[::-1]
    # print(sortedInd)
    # input()
    infoLines = np.array(infoLines)[sortedInd]
    vol = np.array(vol)[sortedInd]

    nPockets = len(infoLines)

    print("Number of pockets found: ",nPockets)
    if nPockets>0:
        print("Largest volume = ",np.amax(vol))
        print("Smallest volume = ",np.amin(vol))

    summaryFile = open(outFolder+structureName+"_infoPockets.txt",'w')
    boolmap = {}
    alreadySeen=set() #to make sure first ranked pockets have priority in writing when atoms in common..
    for ind,line in enumerate(infoLines):
        patoms = [int(l)-1 for l in line.split()]
        # print(patoms)
        # input() 
        summaryFile.write("Pocket %d\n  Vol=%.2f\n"%(ind+1,np.round(vol[ind],2)))
        if(makeMAP):
            #create mapping
            for pindex in patoms:
                if pindex in alreadySeen:
                    pass
                else:
                    boolmap[pindex] = ind+1  
                    alreadySeen.add(pindex)
        else:
            savePQR(patoms,resMap,ind,outFolder)
        if(makeOFF):
            subprocess.run(['mv',runFolder+"cav_tri"+str(sortedInd[ind])+".off",outFolder+'p'+str(ind+1)+".off"])
    

    # print(boolmap)
    # Save structure tringulation
    if(makeOFF):
        subprocess.run(['mv',runFolder+"triangulatedSurf.off",outFolder+structureName+".off"])

    if(makeMAP):
        outMap = open(outFolder+structureName+"_labelled.pqr",'w')
        for ind in range(len(resMap)):
            if ind in boolmap.keys():
                dummyCharge = boolmap[ind]
            else: 
                dummyCharge = 0
            # try:
            #     rChain = resMap[ind]['resChain']
            # except KeyError:
            #     rChain = 'A'
            # outMap.write(str(boolmap[i])+'\n')

            outMap.write("{:<6s}{:>5d} {:<5s}{:>3s} {:>5s}   {:>8.3f} {:>8.3f} {:>8.3f} {:>8.4f} {:>8.4f}\n".format('ATOM',resMap[ind]['atomNumber'],
            resMap[ind]['resAtom'],resMap[ind]['resName'],resMap[ind]['resNum'],
            resMap[ind]['coord'][0],resMap[ind]['coord'][1],resMap[ind]['coord'][2],dummyCharge,resMap[ind]['radius']))
            
    outMap.close()


    end=datetime.now()
    elapsed = end-start
    summaryFile.write("\n---------\nElapsed time: {}".format(elapsed))
    summaryFile.close()
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\nUser exit")
        sys.exit()

