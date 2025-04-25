#### Author: Mikael Trellet (February 2011)

#### Get irmsd and haddock score for ####
#### each complex. You just have to  ####
#### launch the script from your     ####
#### CAPRI_analysis directory (rmsd  ####
#### required in i-RMSD.dat           ####

import os
from os import _exit
from sys      import argv, stderr, stdout
from math         import sqrt
import re

score={} # List of haddock scores (score[1] -> haddock score of complex_1w.pdb)
irmsd={} # List of i-rmsd (irmsd[1] -> i-RMSD of complex_1w.pdb)
lrmsd={} # List of l-rmsd (lrmsd[1] -> l-RMSD of complex_1w.pdb)
fnat={} # List of fnat (fnat[1] -> FNAT of complex_1w.pdb)
vdw={} # List of Van der Waals energy (vdw[1] -> Van der Waals energy of complex_1w.pdb)
elec={} # List of electrostatics (elec[1] -> electrostatic energy of complex_1w.pdb)
air={} # List of AIRs energies (elec[1] -> electrostatic energy of complex_1w.pdb)
clus={} # List of cluster numbers (clus[1] -> cluster of complex_1w.pdb (0 if not clustered))
ilrmsd={} # List of i-l-rmsd (i-l-rmsd[1] -> i-l-RMSD of complex_1w.pdb)
desolv={} # List of desolvation energies (desolv[1] -> desolvation energy of complex_1w.pdb)
match=[]

def getHaddock(file, *refe): # Get haddock scores for further add to specific analysis files (coupled with irmsd, lrmsd and fnat values)
    for l in file.readlines():
        m=re.findall('[0-9]+',l.split()[0])
        num=int(m[len(m)-1])
        if num:
            score[int(num)]=str(l.split()[2])
    file.close()

def getiRMSD(file): # Get irmsd values to further analysis and graphics
    for l in file.readlines():
        if l.split()[0]!='#Interface':
            m=re.findall('[0-9]+',l.split()[0])
            num=int(m[len(m)-1])
            irmsd[num]=str(l.split()[1])
    file.close()

def getlRMSD(file):
    for l in file.readlines():
        if l.split()[0]!='#Ligand':
            m=re.findall('[0-9]+',l.split()[0])
            num=int(m[len(m)-1])
            lrmsd[num]=str(l.split()[1])
    file.close()

def getilRMSD(file):
    #extract ilrmsd from lrmsd in case of FCC clustering
    for l in file.readlines():
        if l.split()[0]!='#Ligand':
            m=re.findall('[0-9]+',l.split()[0])
            num=int(m[len(m)-1])
            ilrmsd[num]=str(l.split()[1])
    file.close()

def getFnat(file, *refe): # Get fnat values to further analysis and graphics
    o=open('../file.nam','r')
    for l in o.readlines():
        if not l.isspace():
            m=re.findall('[0-9]+',l.split()[0])
            num=int(m[len(m)-1])
            match.append(num)
    o.close()
    for line in file.readlines():
        num=match[int(line.split()[1])-2]
        fnat[num]=str(line.split()[2])
    file.close()

def getVdW_Elec_AIR():
    for file in os.listdir(os.getcwd()+"/../"):
        if file!="reference.pdb" and file[-4:]==".pdb":
            o=open('../'+file, 'r')
            m=re.findall('[0-9]+',file)
            num=int(m[len(m)-1])
            for l in o.readlines():
                if len(l.split())>1 and l.split()[1]=="energies:":
                    vdw[int(num)]=str(l.split()[7].strip(","))
                    elec[int(num)]=str(l.split()[8].strip(","))
                    air[int(num)]=str(l.split()[9].strip(","))

def get_desolv():
    for pdb in os.listdir(os.getcwd()+"/../"):
        if pdb!="reference.pdb" and pdb[-4:]==".pdb":
            m=re.findall('[0-9]+', pdb)
            num=int(m[len(m)-1])
            with open('../'+pdb) as o:
                for l in o.readlines():
                    if len(l.split()) > 1 and l.split()[1] == 'Desolvation':
                        desolv[int(num)] = l.split()[3]

def getInter_ligand_RMSD():
    o=open('../analysis/complex_rmsd.disp','r')
    for l in o:
        ilrmsd[match[0]]='0.00000'
        if l.split()[0]=='1':
            num=match[int(l.split()[1])-1]
            ilrmsd[num]=str(l.split()[2])
    o.close()


def getCluster():
    clusteritos_found = []
    for file in os.listdir(os.getcwd()+"/../"):
        if file.startswith('file.nam_clust'):
            s = file
            clusteritos_found.append(s.replace('file.nam_clust','').split('_')[0].split('.')[0])
    clusteritos_found = list(set(clusteritos_found))
    for clusterito in clusteritos_found:
        file = 'file.nam_clust' + clusterito
        o = open('../' + file, 'r')
        #m = re.findall('[0-9]+', file)
        #num = int(m[len(m)-1])
        for l in o.readlines():
            f = re.findall('[0-9]+', l.split()[0])
            struc = int(f[len(f)-1])
            clus[struc] = clusterito


def main(argv):
    try:
        file=open('../file.list-analysis')
    except:
        print("file.list isn't readable, check if it's in your parent directory...")
        return 1
    getHaddock(file)

    try:
        file=open('i-RMSD.dat')
    except:
        print("i-RMSD.dat isn't readable, check if it's in your current directory...")
        return 1
    getiRMSD(file)

    try:
        file=open('l-RMSD.dat')
    except:
        print("l-RMSD.dat isn't readable, check if it's in your current directory...")
        return 1
    getlRMSD(file)

    try:
        file=open('file.nam_fnat')
    except:
        print("file.nam_fnat isn't readable, check if it's in your current directory...")
    getFnat(file)

    getVdW_Elec_AIR()
    get_desolv()
    if os.path.exists("../analysis/complex_rmsd.disp"):
        getInter_ligand_RMSD()
    else:
        print 'FCC clustering performed, using l-RMSD instead of i-l-RMSD'
        try:
            file=open('l-RMSD.dat')
        except:
            print("l-RMSD.dat isn't readable, check if it's in your current directory...")
            return 1
        getilRMSD(file)
    getCluster()
    result=''

    for i in score.iterkeys():
        if i not in clus.keys():
            result=result+(str(i)+" "+str(score[i])+" "+str(irmsd[i])+" "+str(lrmsd[i])+" "+str(fnat[i])+" "+str(vdw[i])+" "+str(elec[i])+" "+str(air[i])+" "+str(desolv[i])+" "+str(ilrmsd[i])+" 0"+"\n")
        else:
            result=result+(str(i)+" "+str(score[i])+" "+str(irmsd[i])+" "+str(lrmsd[i])+" "+str(fnat[i])+" "+str(vdw[i])+" "+str(elec[i])+" "+str(air[i])+" "+str(desolv[i])+" "+str(ilrmsd[i])+" "+str(clus[i])+"\n")

    result_file=open('complex_HS_irmsd_lrmsd_fnat.list','w')
    result_file.write(result)
    result_file.close()

    return 0

end = main(argv)
_exit(end)
