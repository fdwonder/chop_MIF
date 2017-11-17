######################################################################################################################################
# Automatic Zmatrix script generator for BOSS and MCPRO 
# Author: Israel Cabeza de Vaca Lopez   
# Date: March 2015  
# Place: Prof. Jorgensen Lab, Chemistry department at YALE University
# email: israel.cabezadevaca@yale.edu  or israel.cabeza@bsc.es
# 
# Script based on the README notes written by Dr. Julian Tirado-Rives 
#
# Usage:
#   Edit the input parameters written at the beginning of the file for your system. (Variable names are so clear  :) )
# Outputs:
#
#   It generates a folder called finalZmatrices with the final zmatrix files (all, cap, cap + conrot for flexible backbones) 
#
######################################################################################################################################



 
### ********************** INPUT PARAMETERS *********************************

# List of Inputs to perform the calculation

complexPDB = '1wbv.pdb'
ligandZmatOrig = 'benzene.z'
ligandResidueName = 'LI3'

ligandLstToRemoveFromPDB = [] #recommend you start this empty
residueToBeTheCenterOfTheChopping = 'LI3'   # Normally the ligand
setCapOriginAtom = 'LI3'  # LIGAND name ex. 'LI3'   
setCutOffSize = '18.0'
HipLstOfResidues = []  # resnumber, Chain ex. ['77a','56b'] #optional
HieLstOfResidues = []  # resnumber, Chain 

titleOfYourSystem = 'HIV-RT' # optional

fixBackBoneSelection = [] #['4 67 70 74 110 144 152 153'] # If you have a lot of residues split the selection in different lines
fixShortChainCri = '5'

### ********************** END OF INPUT PARAMETERS *********************************




import os
import sys
from biopandas.pdb import PandasPdb
import pandas as pd
import argparse

def FixShortChain(complexPDB,fixShortChainCri):

    global fixBackBoneSelection

    complexPDBName = complexPDB[:-4]
    choppedPDB = complexPDBName+'.chop.pdb'

    with open(choppedPDB) as f:
        ilines = f.readlines() 

    start_idx = 0
    for num, line in enumerate(ilines):
        resnumlisttmp = []
        if ('TER' in line) and (line.split()[1] != ligandResidueName): #make sure not ligand
            #print(num,line)
            for ele in ilines[start_idx:num]:
                resnumlisttmp.append(ele.split()[4]) # append residue num in one chain to the list
            unique_res_nums = (set(resnumlisttmp))
            if len(unique_res_nums)<=float(fixShortChainCri):
                for res_num in unique_res_nums:
                    fixBackBoneSelection.append(res_num)
            start_idx = num + 1

    #sort the list (not truly necessary)
    fixBackBoneSelection = sorted(fixBackBoneSelection, key=int)

    #turn into correct string format
    fixBackBoneSelection_tmp = fixBackBoneSelection
    fixBackBoneSelection = [' '.join(fixBackBoneSelection_tmp)]            

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
                if int(resnumber)>390 and int(resnumber)<999:
                    ligandResnumber = '  1'
                    pdbLigandFile.write(newLine.replace(resnumber,ligandResnumber))
                elif int(resnumber)>999:
                    ligandResnumber = '   1'
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

    # cmd = '/Applications/Chimera.app/Contents/Resources/bin/chimera --nogui protonate.cmd'
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

def generateScript(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HieLstOfResidues):

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

    tmpScript = tmpScript + 'add center :'+residueToBeTheCenterOfTheChopping+'\n'
    tmpScript = tmpScript + 'set cap origin :'+setCapOriginAtom+'\n'

    tmpScript = tmpScript + 'set cut origin ligand \n'
    tmpScript = tmpScript + 'set cut size '+setCutOffSize+'\n'
    #tmpScript += 'delete cut :gol\n'
    if len(ligandLstToRemoveFromPDB)!=0:
        for ele in ligandLstToRemoveFromPDB:
            tmpScript = tmpScript + 'delete cut :'+ele.lower()+'\n'

    tmpScript = tmpScript + 'fix chains\n'
    tmpScript = tmpScript + 'cap all\n '
    
    if len(HipLstOfResidues)!=0:
        lst = ''
        for ele in HipLstOfResidues:
            lst = lst + ':'+ele.strip()+' '
        tmpScript = tmpScript + 'set hip '+lst+'\n\n'


    if len(HieLstOfResidues)!=0:
        lst = ''
        for ele in HieLstOfResidues:
            lst = lst + ':'+ele.strip()+' '
        tmpScript = tmpScript + 'set hip '+lst+'\n\n'

    tmpScript = tmpScript + lastPart.replace('aaaa',complexPDBName) 

    return tmpScript
            
def getRemoveList(complexPDB, removeList):
    # read in the pdb file
    ppdb = PandasPdb()
    ppdb.read_pdb(complexPDB)

    # create the dataframes
    atom_df = ppdb.df['ATOM']
    hetatm_df = ppdb.df['HETATM']
    others_df = ppdb.df['OTHERS']

    # find the last atom num before ligands, HOH, etc. for use later
    last_atom_num= atom_df['atom_number'].iloc[-1]

    #make atom_df so it contains both atoms and hetatms
    atom_df = pd.concat([hetatm_df, atom_df])
    #sort based on atom_number so it's in order
    atom_df = atom_df.sort_values(by=['atom_number'])
    #reset the indicies of the df
    atom_df = atom_df.reset_index(drop=True)

    # find the index of the last atom before ligands, HOH, SO4, etc.
    last_atom_idx = atom_df.loc[atom_df['atom_number']==last_atom_num].index

    # create df of just these ending hetatms
    end_hetatm_df = atom_df.iloc[last_atom_idx[0]+1:]

    # get names of thse residues that need to be deleted
    res_names = set(list(end_hetatm_df['residue_name']))

    # append them to the list of ligands to be removed
    for res_name in res_names:
        if res_name != ligandResidueName:
            removeList.append(res_name)

    #make sure unique
    removeList = list(set(removeList))
    print('remove list:')
    print(removeList)

    return removeList

def prepareReducedChopped(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HieLstOfResidues):

    print('Preparing Reduced Chopped System')

    # update the list of residues to be removed
    complexPDBName = complexPDB[:-4]
    ligandLstToRemoveFromPDB = getRemoveList(complexPDBName+'_cplx.pdb', ligandLstToRemoveFromPDB)

    chopScript = generateScript(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HieLstOfResidues)

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

    for line in open(complexPDBName+'.var.in'):
        newLine = line
        if 'parameter type ALL *' in line: newLine = line+'$ set conrot\n$ set override domain 1-'+lastResidue+'\n' 
        if '$ set fixed backbone' in line:
            newLine = ''
            for ele in fixBackBoneSelection: 
                newLine = newLine + '$ set fixed backbone '+ele + '\n'
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


def fixPDBprotein(complexPDB,ligandResidueName):

    print('Fixing PDB')

    # remove chain in the ligand to avoid errors

    fileoutName = complexPDB[:-4]+'_fixed.pdb'

    fileout = open(fileoutName,'w')

    for line in open(complexPDB):
        newLine = line
        if ligandResidueName in line:
            newLine = line[:21]+' '+line[22:]
        fileout.write(newLine)

    fileout.close()
    
    return fileoutName

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

    Example usage: python pdb2zmat_chop_mcr.py -p 4WR8_ABC_3TX.pdb -r 3TX -c 18.0

    Or, if you already have an optimized z-matrix:

    Usage: python pdb2zmat_chop_mcr.py -p 4WR8_ABC_3TX.pdb -z 3TX_h.z -r 3TX -c 18.0
    
    REQUIREMENTS:
    BOSS (need to set BOSSdir in bashrc and cshrc)
    MCPRO (need to set MCPROdir in bashrc and cshrc)
    Preferably Anaconda python with following modules
    pandas 
    biopandas
    argparse
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
        setCapOriginAtom = str(args.resname)

    if args.cutoff:
        setCutOffSize = str(args.cutoff)

    #  *********************** CODE STARTS *********************************************

    checkParameters()

    if args.zmat:
        _,resnumber, resnumberLigandOriginal = generateLigandPDB(complexPDB,ligandResidueName)
        ligandZmat = fixDummyAtomNamesDirect(ligandZmatOrig)

    else:
        ligandZmat,resnumber,resnumberLigandOriginal = prepareLigandZmatrix(complexPDB,ligandResidueName,mcproPath,BOSSscriptsPath)

    fixedComplexPDB = fixPDBprotein(complexPDB,ligandResidueName)

    print('MERGE')
    mergeLigandWithPDBProtein(mcproPath,fixedComplexPDB,ligandZmat,resnumberLigandOriginal)

    print('CHOP')
    prepareReducedChopped(fixedComplexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HieLstOfResidues)

    FixShortChain(fixedComplexPDB,fixShortChainCri)

    print('FINAL')
    prepareFinalZmatrixWithPEPz(fixedComplexPDB,titleOfYourSystem,ligandZmat,fixBackBoneSelection)

    createZmatrices(fixedComplexPDB,ligandResidueName)

    relaxProteinLigand(fixedComplexPDB,mcproPath)

    generateFinalStructuresWithCAP(fixedComplexPDB)


    #do some file management
    pdb_id = pdb_id = os.path.splitext(os.path.basename(args.pdb))[0]
    ligand_id = str(args.resname).lower()

    try:
        os.makedirs(pdb_id + '_files')
        os.makedirs(ligand_id + '_files')
        os.makedirs('other_output')
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
    other_files_l = ['sum','log','out','lll.txt','plt.pdb','chopScript.csh','protonate.cmd','tmp.txt']
    for file in other_files_l:
        try:
            os.system('mv ' + file + ' ./other_output/')
        except:
            pass


