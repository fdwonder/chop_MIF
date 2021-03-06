# !/usr/bin/env python
######################################################################################################################################
# Automatic Zmatrix script generator for BOSS and MCPRO 
#
# Authors: Israel Cabeza de Vaca Lopez and Matthew Robinson
#
# Script based on the README notes written by Dr. Julian Tirado-Rives 
#
# Usage: python pdb2zmat.py -p pdb_file -r residue_name -c cutoff_size
#
#           use --help for further details or instructions
#
# Outputs:
#       It generates a folder called finalZmatrices with the final zmatrix files (all, cap, cap + conrot for flexible backbones) 
# 
# Requirements:
#       BOSS
#       MCPRO
#       Chimera
#       Reduce executable (http://kinemage.biochem.duke.edu/software/reduce.php)
#       Propka-3.1 executable (https://github.com/jensengroup/propka-3.1)
#       Anaconda Python 3.6
#       The conda choptools-env that can be installed from here (https://github.com/robimc14/choptools)
#           -This program should be run within that environment
#       
######################################################################################################################################



 
### ********************** INPUT PARAMETERS *********************************

# List of Inputs to perform the calculation
# Do not edit unless you know what you are doing
# Please see the bottom of the file for more detailed instructions

complexPDB = ''#'1wbv.pdb'
ligandZmatOrig = ''#'benzene.z'
ligandResidueName = ''#'LI3'

ligandLstToRemoveFromPDB = [] #recommend you start this empty, so I need to get rid of these to start
residueToBeTheCenterOfTheChopping = '' #'LI3'   # Normally the ligand
setCapOriginAtom = '' #'LI3'  # LIGAND name ex. 'LI3'   
setCutOffSize = '' #'18.0'
HipLstOfResidues = []  # resnumber, Chain ex. ['77a','56b'] #optional
#HieLstOfResidues = []  # resnumber, Chain this the default!!! original code wrong
HidLstOfResidues = []

titleOfYourSystem = '' # optional

fixBackBoneSelection = [] #['4 67 70 74 110 144 152 153'] # If you have a lot of residues split the selection in different lines

### ********************** END OF INPUT PARAMETERS *********************************




import os
import sys
import subprocess
from biopandas.pdb import PandasPdb
import pandas as pd
import argparse

def runPropka(originalPDB,HipLstOfResidues):
    #run propka on the original protein


    #try:
    os.system('propka31 ' + originalPDB)
    #except:
        #print('propka failed')
        #return ''

    propka_output_name = complexPDB[:-4]+'.pka'

    #read propka output
    with open(propka_output_name) as propka_output:
        pka_data = propka_output.readlines()

    #get only the pKa data
    for i in range(len(pka_data)):
        if (pka_data[i][0:7] == 'SUMMARY'):
            pka_data = pka_data[i+1:]
            break

    #read the histidine residues
    propka_hip_list = [] 
    for line in pka_data:
        line = line.rstrip()
        if line.split() == []:
            continue
        res_name = line.split()[0]
        if res_name == 'HIS':
            res_number = line.split()[1] + line.split()[2]
            pka = line.split()[3]
            if (float(pka)>=7.0):
                propka_hip_list.append(res_number)

    #compare the two lists
    for res in propka_hip_list:
        if res not in HipLstOfResidues:
            print("prokpka thinks " + res + " should be HIP, reduce does not")

    for res in HipLstOfResidues:
        if res not in propka_hip_list:
            print("reduce thinks " + res + " should be HIP, propka does not")

def makeHisLists(originalPDB,complexPDB,HipLstOfResidues,HidLstOfResidues):

    pdb_name = complexPDB[:-4]+'_cplx.pdb'

    pdb_out = complexPDB[:-4]+'_cplx_h.pdb'

    # run reduce on the protein
    #subprocess.call("reduce –build " + pdb_name + " > " + pdb_name[:-4]+ "_h.pdb", shell=True)
    os.system("reduce -Build %s > %s" % (pdb_name, pdb_out))
        
    
    #read in the protein with biopandas
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_out)
    atom_df = ppdb.df['ATOM']

    #construct variables for boolean selection
    HID = atom_df['atom_name']=='HD1'
    HIE = atom_df['atom_name']=='HE2'
    HIS = atom_df['residue_name']=='HIS'
    
    #make dataframes
    hid_df = atom_df.loc[HIS & HID]
    hie_df = atom_df.loc[HIS & HIE]
    
    #make booleans to check it has that type of his
    has_hid = hid_df.shape[0] > 0
    has_hie = hie_df.shape[0] > 0
    
    #construct the lists
    hid_list = []
    hie_list = []
    hip_list = []
    
    #first for his
    if has_hid:
        hid_nums = [str(x) for x in list(hid_df['residue_number'])]
        hid_chains = list(hid_df['chain_id'])
        hid_list = [hid_nums[i] + hid_chains[i] for i in range(len(hid_nums))]
    #now for hie
    if has_hie:
        hie_nums = [str(x) for x in list(hie_df['residue_number'])]
        hie_chains = list(hie_df['chain_id'])
        hie_list = [hie_nums[i] + hie_chains[i] for i in range(len(hie_nums))]
    
    #construct the hip list
    if has_hid and has_hie:
        for res in hie_list:
            if res in hid_list:
                hip_list.append(res)

    print('hiplist')
    print(hip_list)
    print('hielist')
    print(hie_list)
    print('hidlist')
    print(hid_list)

    #add to the global lists
    for res in hip_list:
        if (res not in HipLstOfResidues) and (res not in HidLstOfResidues):
            HipLstOfResidues.append(res)

    for res in hid_list:
        if (res not in HipLstOfResidues) and (res not in HidLstOfResidues):
            HidLstOfResidues.append(res)

    #check if propka disagrees with reduce
    runPropka(originalPDB,HipLstOfResidues)

    return HipLstOfResidues, HidLstOfResidues     

def checkParameters():

    if not mcproPath:
        print('Define MCPROdir enviroment variable, please')
        sys.exit()
    else:
        print('Using MCPROdir .....    ',mcproPath)

    if not BOSSPath:
        print('Define BOSSdir enviroment variable, please')
        sys.exit()
    else:
        print('Using BOSSdir .....    ',BOSSPath)

    if not os.path.isfile(complexPDB):
        print('PDB not found ......    ',complexPDB)
        sys.exit()

def generateLigandPDB(complex,ligandName):

    print('Generating Ligand PDB')

    namePDBLigand = ligandName.lower()+'.pdb'
    pdbLigandFile = open(namePDBLigand,'w')

    ligandResnumber = ''

    ligandResnumberOriginal = ''

    for line in open(complex):
        if 'ATOM' in line[:7] or 'HETATM' in line[:7]:
            if ligandName in line[15:23]:
                
                newLine = line[:21]+' '+line[22:]
                resnumber = newLine[12:29].split()[2]
                ligandResnumber = resnumber
                ligandResnumberOriginal = resnumber
                if int(resnumber)<=999:
                    ligandResnumber = '100'
                    pdbLigandFile.write(newLine.replace(resnumber,ligandResnumber))
                elif int(resnumber)>999:
                    ligandResnumber = ' 100'
                    pdbLigandFile.write(newLine.replace(resnumber,ligandResnumber))
                else:
                    pdbLigandFile.write(newLine)

    pdbLigandFile.close()

    #complexPDBfixName = fixComplexResidueNumber(complex,ligandName)

    return namePDBLigand,ligandResnumber,ligandResnumberOriginal#,complexPDBfixName
    

chimeraScript='''open aaaa
addh
write 0 bbbb

'''

def protonateLigandWithChimera(namePDBLigand):

    print('Protonating ligand with Chimera')

    namePDBLigand_protonated = namePDBLigand[:-4]+'_h.pdb'

    chimeraScriptFile = open('protonate.cmd','w')
    chimeraScriptFile.write(chimeraScript.replace('aaaa',namePDBLigand).replace('bbbb',namePDBLigand_protonated))
    chimeraScriptFile.close()

    #cmd = '/Applications/Chimera.app/Contents/Resources/bin/chimera --nogui protonate.cmd'
    cmd = 'chimera --nogui protonate.cmd'

    os.system(cmd)

    return namePDBLigand_protonated 

def fixDummyAtomNames(namePDBLigand_protonated):
    
    print('Fixing DUMMY atoms names')

    filetmp = open('tmp.txt','w')

    ligandZmat = namePDBLigand_protonated[:-4]+'.z'

    for line in open(ligandZmat):
        newLine = line
                
        if len(line.split())==12 and len(line)==77:
            newLine = line[:71]+'    1\n'
        filetmp.write(newLine.replace('   1 DUM','   1 DU1').replace('   2 DUM','   2 DU2').replace('   3 O00  800  800','   3 DU3   -1    0'))

    filetmp.close()

    cmd = 'cp tmp.txt '+ligandZmat
    
    print(cmd)

    os.system(cmd)

    return ligandZmat

def fixDummyAtomNamesDirect(ligandZmatOrig):

    print('Fixing DUMMY atoms names')

    filetmp = open('tmp.txt', 'w')

    ligandZmat = ligandZmatOrig[:-2] + '_fixed.z'

    for line in open(ligandZmatOrig):
        newLine = line

        if len(line.split()) == 12 and len(line) == 77:
            newLine = line[:71] + '    1\n'
        filetmp.write(
            newLine.replace('   1 DUM', '   1 DU1').replace('   2 DUM', '   2 DU2').replace('   3 O00  800  800',
                                                                                            '   3 DU3   -1    0'))

    filetmp.close()

    cmd = 'cp tmp.txt ' + ligandZmat

    print(cmd)

    os.system(cmd)

    return ligandZmat

def optimizeStructure(ligandZmat,mcproPath):

    print('Optimizing Structure')

    cmd = mcproPath+'/scripts/xOPT '+ligandZmat[:-2]

    print(cmd)

    os.system(cmd)

def prepareLigandZmatrix(complex,ligandName,mcproPath,BOSSscriptsPath):

    print('Preparing Ligand Z matrix')

    namePDBLigand, ligandResnumber, ligandResnumberOriginal = generateLigandPDB(complex,ligandName)

    namePDBLigand_protonated = protonateLigandWithChimera(namePDBLigand)

    cmd = BOSSscriptsPath+'/xPDBMCP '+namePDBLigand_protonated[:-4]

    os.system(cmd)

    ligandZmat = fixDummyAtomNames(namePDBLigand_protonated)
    
    optimizeStructure(ligandZmat, mcproPath)

    return ligandZmat,ligandResnumber,ligandResnumberOriginal#,complexPDBfixName

def mergeLigandWithPDBProtein(mcproPath,complexPDB,ligandResidueName,resnumber):
    #clu -t:s=5001 2be2.pdb -r r22_h.z -n 2be2_cplx.pdb

    print('Merging ligand with PDB Protein')
    
    clu = mcproPath+'/miscexec/clu '

    outputPDB = complexPDB[:-4]+'_cplx.pdb'
    
    if os.path.isfile(outputPDB): os.remove(outputPDB)

    cmd = clu + ' -t:s='+resnumber+' '+complexPDB+' -r '+ligandZmat+' -n '+outputPDB

    print(cmd)

    os.system(cmd)

def generateScript(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HidLstOfResidues):

    lastPart = '''set variable origin ligand
set variable size 10  
write pdb aaaa.chop.pdb
write pepz all aaaa.chop.all.in
write pepz variable aaaa.chop.var.in
write translation aaaa.chop.tt
EOF

'''

    complexPDBName = complexPDB[:-4]

    tmpScript = '$MCPROdir/miscexec/chop -t -i '+complexPDBName+'_cplx.pdb << EOF\n'

    if len(ligandLstToRemoveFromPDB)!=0:
        for ele in ligandLstToRemoveFromPDB:
            tmpScript = tmpScript + 'delete ligand :'+ele+'\n'

    if len(HipLstOfResidues)!=0:
        lst = ''
        for ele in HipLstOfResidues:
            lst = lst + ':'+ele.strip()+' '
        tmpScript = tmpScript + 'set hip '+lst+'\n\n'


    if len(HidLstOfResidues)!=0:
        lst = ''
        for ele in HidLstOfResidues:
            lst = lst + ':'+ele.strip()+' '
        tmpScript = tmpScript + 'set hid '+lst+'\n\n'

    tmpScript = tmpScript + 'add center :'+residueToBeTheCenterOfTheChopping+'\n'
    tmpScript = tmpScript + 'set cap origin :'+setCapOriginAtom+'\n'

    tmpScript = tmpScript + 'set cut origin ligand \n'
    tmpScript = tmpScript + 'set cut size '+setCutOffSize+'\n'
    #tmpScript += 'delete cut :gol\n'
    if len(ligandLstToRemoveFromPDB)!=0:
        for ele in ligandLstToRemoveFromPDB:
            tmpScript = tmpScript + 'delete cut :'+ele.lower()+'\n'

    tmpScript = tmpScript + 'set minchain 5\n'
    tmpScript = tmpScript + 'fix chains\n'
    tmpScript = tmpScript + 'cap all\n '
    
    # if len(HipLstOfResidues)!=0:
    #     lst = ''
    #     for ele in HipLstOfResidues:
    #         lst = lst + ':'+ele.strip()+' '
    #     tmpScript = tmpScript + 'set hip '+lst+'\n\n'


    # if len(HidLstOfResidues)!=0:
    #     lst = ''
    #     for ele in HidLstOfResidues:
    #         lst = lst + ':'+ele.strip()+' '
    #     tmpScript = tmpScript + 'set hid '+lst+'\n\n'

    tmpScript = tmpScript + lastPart.replace('aaaa',complexPDBName) 

    return tmpScript

def prepareReducedChopped(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HidLstOfResidues):

    print('Preparing Reduced Chopped System')

    # update the list of residues to be removed
    complexPDBName = complexPDB[:-4]

    chopScript = generateScript(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HidLstOfResidues)

    fileChopScript = open('chopScript.csh','w')
    fileChopScript.write(chopScript)
    fileChopScript.close()

    os.system('csh chopScript.csh')

def getNumberOfTheLastResidueOfTheChoppedSystem(choppedPDBName):
    
    resNumber = ''

    for line in open(choppedPDBName):
        if 'ATOM' in line[:6] or 'HETATM' in line[:6]:
            resNumber = line.split()[4]

    return resNumber    


def prepareFinalZmatrixWithPEPz(complexPDB,titleOfYourSystem,ligandZmat,fixBackBoneSelection):

    complexPDBName = complexPDB[:-4]

    tmpfile = open(complexPDBName+'.all.in','w')

    for line in open(complexPDBName+'.chop.all.in'):
        newLine = line
        if '[ADD YOUR TITLE HERE]' in line: newLine = line.replace('[ADD YOUR TITLE HERE]',titleOfYourSystem)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: newLine = line.replace('[WRITE NAME OF YOUR solute z-matrix FILE]',ligandZmat)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: newLine = line.replace('[NAME OF THE z-matrix TO BE WRITTEN]',complexPDBName+'all.z')

        tmpfile.write(newLine)

    tmpfile.close()
        
    tmpfile = open(complexPDBName+'.var.in','w')

    for line in open(complexPDBName+'.chop.var.in'):
        newLine = line
        if '[ADD YOUR TITLE HERE]' in line: newLine = line.replace('[ADD YOUR TITLE HERE]',titleOfYourSystem)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: newLine = line.replace('[WRITE NAME OF YOUR solute z-matrix FILE]',ligandZmat)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: newLine = line.replace('[NAME OF THE z-matrix TO BE WRITTEN]',complexPDBName+'var.z')

        tmpfile.write(newLine)

    tmpfile.close()

    # create conrot

    lastResidue = getNumberOfTheLastResidueOfTheChoppedSystem(complexPDBName+'.chop.pdb')
    
    tmpfile = open(complexPDBName+'.var_conrot.in','w')

    #get the residues that need to be fixed for conrot
    res_strings = []
    for line in open(complexPDBName+'.var.in'):
        #get the initial strings of the fixed res
        if 'set fixed all' in line:
            res_strings = res_strings + line.split()[4:]

    #get all the numbers
    res_nums = []
    for res_str in res_strings:
        if '-' in res_str:
            flanks = res_str.split('-')
            res_nums = res_nums + list(range(int(flanks[0]),int(flanks[1])+1))
        else:
            res_nums = res_nums + [int(res_str)]

    #make the list of residues to fix the backbone of 
    fix_list = []
    for i in range(1,len(res_nums)):
        num = res_nums[i]
        prev_num = res_nums[i-1]
        diff = num - prev_num
        if (diff <= 5) and (diff != 1):
            fix_list = fix_list + list(range(prev_num+1,num))

    #put list in correct str format        
    fix_list_tmp = [str(x) for x in fix_list]
    fix_list = [' '.join(fix_list_tmp)]        


    for line in open(complexPDBName+'.var.in'):
        newLine = line
        if 'parameter type ALL *' in line: newLine = line+'$ set conrot\n$ set override domain 1-'+lastResidue+'\n' 
        if '$ set fixed backbone' in line:
            newLine = ''
            #for ele in fixBackBoneSelection: 
            for ele in fix_list:
                newLine = newLine + '$ set fixed backbone '+str(ele) + '\n'
        if '$ write zmatrix '+complexPDBName+'var.z' in line: newLine = line.replace(complexPDBName+'var.z',complexPDBName+'varcon.z')

        tmpfile.write(newLine)

    tmpfile.close()

def fixZmatrix(matrixZ,ligandResidueName):

    print('Fixing Z matrix TERZ .... ',matrixZ)

    tmpfile = open('tmp.txt','w')

    be2allFile = [ele for ele in open(matrixZ)]
    
    for iter in range(len(be2allFile)):
        newLine = be2allFile[iter]
        if 'TERZ' in newLine:
            if ligandResidueName in be2allFile[iter+1] or 'CAP' in be2allFile[iter+1]:
                newLine = be2allFile[iter]
            else:
                continue

        tmpfile.write(newLine)

    tmpfile.close()
                

    os.system('cp '+matrixZ+' lll.txt') 
    os.system('cp tmp.txt '+matrixZ)    
    
def createZmatrices(complexPDB,ligandResidueName):

    MCPROscriptsPath = os.path.join(mcproPath,'scripts')

    complexPDBName = complexPDB[:-4]

    os.system(MCPROscriptsPath+'/xPEPZ '+complexPDBName+'.all')
    os.system(MCPROscriptsPath+'/xPEPZ '+complexPDBName+'.var')
    os.system(MCPROscriptsPath+'/xPEPZ '+complexPDBName+'.var_conrot')

    fixZmatrix(complexPDBName+'all.z',ligandResidueName)        
    fixZmatrix(complexPDBName+'var.z',ligandResidueName)        
    fixZmatrix(complexPDBName+'varcon.z',ligandResidueName)     

def relaxProteinLigand(complexPDB,mcproPath):
    
    complexPDBName = complexPDB[:-4]

    os.system('mkdir CG9;cp '+complexPDBName+'all.z CG9;cd CG9;'+mcproPath+'/scripts/xCGDD9 '+complexPDBName+'all')
            

def getOptimizedZmatrix():

    # By default the file is called optsum 

    ZmatrixOptimized = []

    for line in open('CG9/optsum'):

        if 'Geometry Variations follow' in line: break

        ZmatrixOptimized.append(line)

    return ZmatrixOptimized

def replaceOptimizedCoordinatesInZmatrixFile(varFileName,capFileName,zmatrixOptimized):

    fileout = open(capFileName,'w')

    # get bottom part of the var file

    read = False
    bottomDataInVarFile = []
    
    for line in open(varFileName):
        if 'Geometry Variations follow' in line: read = True
        if read: bottomDataInVarFile.append(line)

    for ele in zmatrixOptimized: fileout.write(ele)
    for ele in bottomDataInVarFile: fileout.write(ele)

    
def generateFinalStructuresWithCAP(complexPDB):

    print('Generating final structures')

    complexPDBName = complexPDB[:-4]

    zmatrixOptimized = getOptimizedZmatrix()


    replaceOptimizedCoordinatesInZmatrixFile(complexPDBName+'var.z',complexPDBName+'cap.z',zmatrixOptimized)
    replaceOptimizedCoordinatesInZmatrixFile(complexPDBName+'varcon.z',complexPDBName+'capcon.z',zmatrixOptimized)
    
    try:
        os.mkdir('finalZmatrices')
    except:
        os.system('rm -r finalZmatrices')
        os.mkdir('finalZmatrices')

    os.system('cp CG9/optsum finalZmatrices/'+complexPDBName+'all.z')   
    os.system('cp '+complexPDBName+'cap.z '+ complexPDBName+'capcon.z finalZmatrices')


def removeSolvent(complexPDB):

    print('Removing Solvent')

    fileoutName = complexPDB[:-4]+'_no_solvent.pdb'

    ppdb = PandasPdb().read_pdb(complexPDB)
    atom_df = ppdb.df['ATOM']

    # could make this a for loop with the delete lists
    atom_df = atom_df[atom_df['residue_name'] != 'HOH']
    atom_df = atom_df[atom_df['residue_name'] != 'SO4']

    ppdb.to_pdb(path=fileoutName, 
            records=None, 
            gz=False, 
            append_newline=True)
    
    return fileoutName

def fixPDBprotein(complexPDB,ligandResidueName):

    print('Fixing PDB')

    # remove chain in the ligand to avoid errors

    fileoutName = complexPDB[:-4]+'_fixed.pdb'

    fileout = open(fileoutName,'w')

    for line in open(complexPDB[:-4]+'_no_solvent.pdb'):
        newLine = line
        if ligandResidueName in line:
            newLine = line[:21]+' '+line[22:]
        fileout.write(newLine)

    fileout.close()
    
    return fileoutName

def manageFiles():
    pdb_id = pdb_id = os.path.splitext(os.path.basename(args.pdb))[0]
    ligand_id = str(args.resname).lower()

    try:
        os.makedirs(pdb_id + '_files')
        os.makedirs(ligand_id + '_files')
    except:
        print('folders already made')
        pass

    #move the pdb files
    pdb_output_files = [filename for filename in os.listdir('.') if filename.startswith(pdb_id)]
    #remove the folder name
    if (pdb_id + '_files') in pdb_output_files:
        pdb_output_files.remove(pdb_id + '_files')
    #mv these files to the output folder
    for file in pdb_output_files:
        try:
            os.system('mv ' + file + ' ./' + pdb_id + '_files/')   
        except:
            pass

    #now same thing for the ligand
    ligand_output_files = [filename for filename in os.listdir('.') if filename.startswith(ligand_id)]
    #remove the folder name
    if (ligand_id + '_files') in ligand_output_files:
        ligand_output_files.remove(ligand_id + '_files')
    #mv these files to the output folder
    for file in ligand_output_files:
        try:
            os.system('mv ' + file + ' ./' + ligand_id + '_files/')   
        except:
            pass

    #take care of the other files
    try:
        os.makedirs(pdb_id + 'other_output')
    except:
        pass
    other_files_l = ['sum','log','out','lll.txt','plt.pdb','chopScript.csh','protonate.cmd','tmp.txt']
    for file in other_files_l:
        try:
            os.system('mv ' + file + ' ./' + pdb_id + 'other_output/')
        except:
            pass

if __name__ == "__main__":

    # get paths of needed dirs
    mcproPath = os.environ.get('MCPROdir')
    BOSSPath = os.environ.get('BOSSdir')

    BOSSscriptsPath = os.path.join(BOSSPath,'scripts')

    # create parser object
    parser = argparse.ArgumentParser(
        prog='pdb2zmatConrotGenerator_mcr.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
    Automatic Zmatrix script generator for BOSS and MCPRO 

    @author: Israel Cabeza de Vaca Lopez, israel.cabezadevaca@yale.edu   
    @author: Matthew Robinson, matthew.robinson@yale.edu
    @author: Yue Qian, yue.qian@yale.edu
    @author: William L. Jorgensen Lab 

    Example usage: python pdb2zmat.py -p 4WR8_ABC_3TX.pdb -r 3TX -c 18.0

    Or, if you already have an optimized z-matrix:

    Usage: python pdb2zmat.py -p 4WR8_ABC_3TX.pdb -z 3TX_h.z -r 3TX -c 18.0
    
    REQUIREMENTS:
    BOSS (need to set BOSSdir in bashrc and cshrc)
    MCPRO (need to set MCPROdir in bashrc and cshrc)
    Chimera
    reduce executable (from Richardson lab)
    propka-3.1 executable (from Jensen lab)
    Preferably Anaconda python 3.6 with following modules:
    pandas 
    biopandas
    """
    )

    #defining arguments for parser object
    parser.add_argument(
        "-p", "--pdb", help="name of the PDB file for the complex - i.e., [complex_name].pdb")
    parser.add_argument(
        "-z", "--zmat", help="name of the zmat of the ligand with .z file descriptor. \
        only need this if want to import the optimized ligand zmat yourself.")
    parser.add_argument(
        "-r", "--resname", help="Residue name of the ligand from PDB FILE", type=str)
    parser.add_argument(
        "-c", "--cutoff", help="Size of the cutoff for chop to cut", type=str)

    #parse the arguments from standard input
    args = parser.parse_args()

    if args.zmat:
        ligandZmatOrig = str(args.zmat)

    if args.pdb:
        complexPDB = str(args.pdb)

    if args.resname:
        ligandResidueName = str(args.resname)
        residueToBeTheCenterOfTheChopping = str(args.resname) 
        setCapOriginAtom = 'c01'#str(args.resname)

    if args.cutoff:
        setCutOffSize = str(args.cutoff)

    #  *********************** CODE STARTS *********************************************

    checkParameters()

    if args.zmat:
        _,resnumber, resnumberLigandOriginal = generateLigandPDB(complexPDB,ligandResidueName)
        ligandZmat = fixDummyAtomNamesDirect(ligandZmatOrig)

    else:
        ligandZmat,resnumber,resnumberLigandOriginal = prepareLigandZmatrix(complexPDB,ligandResidueName,mcproPath,BOSSscriptsPath)


    #remove solvent
    no_solvent_PDB = removeSolvent(complexPDB)

    #change Chain ID of ligand
    fixedComplexPDB = fixPDBprotein(complexPDB,ligandResidueName)


    print('MERGE')
    mergeLigandWithPDBProtein(mcproPath,fixedComplexPDB,ligandZmat,resnumberLigandOriginal)

    #get histdine lists before chopping
    HipLstOfResidues, HidLstOfResidues = makeHisLists(args.pdb, fixedComplexPDB,HipLstOfResidues,HidLstOfResidues)
    print(HipLstOfResidues)
    print(HidLstOfResidues)

    print('CHOP')
    prepareReducedChopped(fixedComplexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HidLstOfResidues)

    print('FINAL')
    prepareFinalZmatrixWithPEPz(fixedComplexPDB,titleOfYourSystem,ligandZmat,fixBackBoneSelection)

    createZmatrices(fixedComplexPDB,ligandResidueName)

    relaxProteinLigand(fixedComplexPDB,mcproPath)

    generateFinalStructuresWithCAP(fixedComplexPDB)

    #do some file management
    manageFiles()

    