# ********************** INPUT PARAMETERS *********************************

# List of Inputs to perform the calculation

complexPDB = 'full.lys.plus8.OPT-orig.pdb'
ligandZmatOrig = 'benzene.z'
ligandResidueName = 'LIG'

ligandLstToRemoveFromPDB = []
residueToBeTheCenterOfTheChopping = 'LIG'   # Normally the ligand
setCapOriginAtom = 'C00'  # the first one before DUMMY atoms in the ligand (normally)
setCutOffSize = '20.0'
HipLstOfResidues = []  # resnumber, Chain  #optional
HidLstOfResidues = []  # resnumber, Chain

titleOfYourSystem = 'LYS-BEN-20A' # optional
ResidueChangeSelection = ['PROA1','PROB1','PROC1']
ResidueFinalSelection = ['PRN', 'PRN', 'PRN'] # ResFiSel has to be of same length as ResChSel

fixBackBoneSelection = [' '] #['4 67 70 74 110 144 152 153'] # If you have a lot of residues split the selection in different lines
fixShortChainCri = '16'

# ********************** END OF INPUT PARAMETERS *********************************
resnumberlist=[]
reslinelist=[]



import os
import sys
import pickle


def locateRes(complexPDB, ResidueChangeSelection,ResidueFinalSelection):
    if len(ResidueChangeSelection)!=0 and len(ResidueFinalSelection)!=0:
        print "Using input files generated from chop, generating output files for pepz"
    else:
        print "Residue change/final list not specified, exiting"
        sys.exit()
    complexPDBName = complexPDB[:-4]+'_cplx.pdb'

    # locate coordinates from complexPDBNamd_fixed_cplx.pdb
    x=[];y=[];z=[]
    for iter in range(len(ResidueChangeSelection)):
        changeres = ResidueChangeSelection[iter]
#        print "now searching", changeres

        for line in open(complexPDBName):
            # match residue name, chain, and residue number and CA only
            if changeres[:3] in line[17:20] and changeres[3:4] in line[21:22] \
            and changeres[4:] == line.split()[5] and 'CA' == line.split()[2]:
                x.append(line.split()[6])
                y.append(line.split()[7])
                z.append(line.split()[8])
    # locate residue number and line in complexPDB_fixed.chop.pdb
    choppedPDB = complexPDB[:-4]+'.chop.pdb'

    resnumberlist = []
    reslinelist =[]
    with open(choppedPDB) as f:
        for num,line in enumerate(f,1):
            newLine=line
            if 'CA' in line[:16]:
                for iter in range(len(ResidueChangeSelection)):
                    if x[iter] == line.split()[5] and y[iter] == line.split()[6] and z[iter] == line.split()[7]:
                        # find the res list from modified chop.pdb file
                        resnumberlist.append(line.split()[4])
                        reslinelist.append(num)

        # print resnumberlist,reslinelist

    # modify chopped pdb file first and name as residue replaced chopped file
    tmpfile = open('tmp.txt', 'w')
    with open(choppedPDB) as f:
        for line in f:
            newLine = line
            for iter in range(len(ResidueChangeSelection)):
                if len(line.split())>=5 and line.split()[4] == resnumberlist[iter]:
                    newLine = line[:17] + ResidueFinalSelection[iter] + line[20:]

            tmpfile.write(newLine)
    tmpfile.close()
    os.system('cp tmp.txt ' + complexPDB[:-4]+'.chop.pdb')
    return resnumberlist,reslinelist


def ReplaceRes(complexPDB,resnumberlist,reslinelist,ResidueFinalSelection):

    complexPDBName = complexPDB[:-4]
    choppedPDB = complexPDBName+'.chop.pdb'
    # count the occurence of TER in chopped pdb file
    flag=[]
    with open(choppedPDB) as f:
        ilines=f.readlines() # f.readlines command is important to have access to line numbers
        for iter in range(len(ResidueChangeSelection)):
            n = 0
            if iter == 0:
                lines =ilines[:reslinelist[0]-1]
            else:
                lines = ilines[reslinelist[iter-1]-1:reslinelist[iter]-1]
                # print 'check last line', ilines[reslinelist[iter]-1]
            for line in lines:
                if 'TER' in line[:3]:
                    n = n + 1
            # append the occurence of TER inbetween each two designated residues to flag
            flag.append(n)

    # replace the residue in all.in , var.in and var_conrot.in files.
    #tmparray = [] #!!! Do not define temparray as list, since you will just need a string
    init = 0; final =0
    with open(complexPDBName+'.chop.all.in') as f:
        for num, line in enumerate(f,1):
            if 'sequence' in line:
                init = num
            if 'center' in line:
                final = num
    # A separate loop is required to read in lines with line number
    with open(complexPDBName+'.chop.all.in') as f:
        ilines = f.readlines()
        tmparray =ilines[init:final-1]
    for i in range(len(ResidueChangeSelection)):
        if i == 0:
            n = int(resnumberlist[i])+flag[i]-1
        else:
            n = n + int(resnumberlist[i])-int(resnumberlist[i-1])+flag[i]

        NumRes=0
        for j in range(len(tmparray)):
            # print 'length of each tmparray ele in j', j, tmparray[j],len(tmparray[j].split())
            oldNumRes = NumRes
            NumRes=NumRes+len(tmparray[j].split())
            if n<NumRes and n >=oldNumRes:
                #str = ''.join(str(ele) for ele in tmparray[j])
                str = tmparray[j] # save the designated line to a string
                A = str.split() # save the splitted string to new list A
                A[n-oldNumRes] = ResidueFinalSelection[i] # assign new residue to the element of the list
                B = ' '.join(A) # convert the list to string
                C = B+'\n'
                tmparray[j] = C # replace the original line with the new line
                # print tmparray

    # Replace corresponding residue in all.in and var.in file
    with open(complexPDBName+'.chop.all.in') as f:
        ilines = f.readlines()
        topprt=ilines[:init]
        botprt=ilines[final-1:]
    with open('tmp.txt', 'w') as tmpfile:
    #      pickle.dump(topprt,tmpfile)
    # for item in topprt:
    #     tmpfile.write('' % item)
        tmpfile.write(''.join(topprt))
        tmpfile.write(''.join(tmparray))
        tmpfile.write(''.join(botprt))
    tmpfile.close()
    os.system('cp tmp.txt ' + complexPDBName +'.chop.all.in')

    with open(complexPDBName+'.chop.var.in') as f:
        ilines = f.readlines()
        topprt=ilines[:init-1]
        botprt=ilines[final-2:]
    with open('tmp.txt', 'w') as tmpfile:
    #      pickle.dump(topprt,tmpfile)
    # for item in topprt:
    #     tmpfile.write('' % item)
        tmpfile.write(''.join(topprt))
        tmpfile.write(''.join(tmparray))
        tmpfile.write(''.join(botprt))
    tmpfile.close()

    os.system('cp tmp.txt ' + complexPDBName +'.chop.var.in')

def FixShortChain(complexPDB,fixShortChainCri):

    complexPDBName = complexPDB[:-4]
    choppedPDB = complexPDBName+'.chop.pdb'
    # count the occurence of TER in chopped pdb file
    flag1=0
    flag2=0
    resnamelisttmp=[]
    resnumlisttmp =[]
    # choppedPDB = 'MIF-180.cm5_fixed.chop.pdb'
    with open(choppedPDB) as f:
        ilines = f.readlines()  # f.readlines command is important to have access to line numbers
    with open(choppedPDB) as f:
        for num, line in enumerate(f, 1):
            if 'TER' in line:
                flag1 = flag1 +1 # occurance of TER
                for ele in ilines[flag2:num-1]:
                    resnamelisttmp.append(ele.split()[3]) # append residue names in one chain to the list
                A = resnamelisttmp
                B = len(set(A)) # count the number of residues on the chain
                if B<=float(fixShortChainCri):
                    for ele in ilines[flag2:num-1]:
                        resnumlisttmp.append(ele.split()[4])  # append residue names in one chain to the list
                    C = set(resnumlisttmp)
                    CC = ' '.join(C)
                    fixBackBoneSelection.append(CC)
                    # fixBackBoneSelection = ' '.join(fixBackBoneSelection)
                flag2 = num


def checkParameters():

	if not mcproPath:
		print 'Define MCPROdir enviroment variable, please'
		sys.exit()
	else:
		print 'Using MCPROdir .....    ',mcproPath

	if not BOSSPath:
		print 'Define BOSSdir enviroment variable, please'
		sys.exit()
	else:
		print 'Using BOSSdir .....    ',BOSSPath

	if not os.path.isfile(complexPDB):
		print 'PDB not found ......    ',complexPDB
		sys.exit()

def generateLigandPDB(complexPDB, ligandResidueName):

    print 'Generating Ligand PDB'

    namePDBLigand = ligandResidueName.lower() + '.pdb'
    pdbLigandFile = open(namePDBLigand, 'w')

    ligandResnumber = ''

    ligandResnumberOriginal = ''

    for line in open(complexPDB):
        if 'ATOM' in line[:7] or 'HETATM' in line[:7]:
            if ligandResidueName in line[15:23]:

                newLine = line[:21] + ' ' + line[22:]
                resnumber = newLine[12:29].split()[2]
                ligandResnumber = resnumber
                ligandResnumberOriginal = resnumber
                if int(resnumber) > 390 and int(resnumber) < 999:
                    ligandResnumber = '  1'
                    pdbLigandFile.write(newLine.replace(resnumber, ligandResnumber))
                elif int(resnumber) > 999:
                    ligandResnumber = '   1'
                    pdbLigandFile.write(newLine.replace(resnumber, ligandResnumber))
                else:
                    pdbLigandFile.write(newLine)

    pdbLigandFile.close()

    # complexPDBfixName = fixComplexResidueNumber(complex,ligandName)

    return namePDBLigand, ligandResnumber, ligandResnumberOriginal  # ,complexPDBfixName

chimeraScript = '''open aaaa
addh
write 0 bbbb

'''

def protonateLigandWithChimera(namePDBLigand):

    print 'Protonating ligand with Chimera'

    namePDBLigand_protonated = namePDBLigand[:-4] + '_h.pdb'

    chimeraScriptFile = open('protonate.cmd', 'w')
    chimeraScriptFile.write(chimeraScript.replace('aaaa', namePDBLigand).replace('bbbb', namePDBLigand_protonated))
    chimeraScriptFile.close()

    cmd = 'chimera --nogui protonate.cmd'

    os.system(cmd)

    return namePDBLigand_protonated


def fixDummyAtomNames(namePDBLigand_protonated):

    print 'Fixing DUMMY atoms names'

    filetmp = open('tmp.txt', 'w')

    ligandZmat = namePDBLigand_protonated[:-4] + '.z'

    for line in open(ligandZmat):
        newLine = line

        if len(line.split()) == 12 and len(line) == 77:
            newLine = line[:71] + '    1\n'
        filetmp.write(
            newLine.replace('   1 DUM', '   1 DU1').replace('   2 DUM', '   2 DU2').replace('   3 O00  800  800',
                                                                                            '   3 DU3   -1    0'))

    filetmp.close()

    cmd = 'cp tmp.txt ' + ligandZmat

    print cmd

    os.system(cmd)

    return ligandZmat

# not part of the prepareLigand Zmatrix, but fix dummy atom names directly from optimized zmat of BOSS
# rename the starting two DUM to DU1 and DU2, duplicate atom names cause problem
def fixDummyAtomNamesDirect(ligandZmatOrig):

    print 'Fixing DUMMY atoms names'

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

    print cmd

    os.system(cmd)

    return ligandZmat


def optimizeStructure(ligandZmat, mcproPath):

    print 'Optimizing Structure'

    cmd = mcproPath + '/scripts/xOPT ' + ligandZmat[:-2]

    print cmd

    os.system(cmd)

def prepareLigandZmatrix(complexPDB, ligandResidueName, mcproPath, BOSSscriptsPath):

    print 'Preparing Ligand Z matrix'

    namePDBLigand, ligandResnumber, ligandResnumberOriginal = generateLigandPDB(complexPDB, ligandResidueName)

    namePDBLigand_protonated = protonateLigandWithChimera(namePDBLigand)

    cmd = BOSSscriptsPath + '/xPDBMCP ' + namePDBLigand_protonated[:-4]

    os.system(cmd)

    ligandZmat = fixDummyAtomNames(namePDBLigand_protonated)

    optimizeStructure(ligandZmat, mcproPath)

    return ligandZmat, ligandResnumber, ligandResnumberOriginal  # ,complexPDBfixName

def fixPDBprotein(complexPDB, ligandResidueName):

    print 'Fixing PDB'

    # remove chain in the ligand to avoid errors

    fileoutName = complexPDB[:-4] + '_fixed.pdb'

    fileout = open(fileoutName, 'w')

    for line in open(complexPDB):
        newLine = line
        if ligandResidueName in line:
            newLine = line[:21] + ' ' + line[22:]
        fileout.write(newLine)

    fileout.close()

    return fileoutName

def mergeLigandWithPDBProtein(mcproPath, complexPDB, ligandResidueName, resnumber):
    # clu -t:s=5001 2be2.pdb -r r22_h.z -n 2be2_cplx.pdb

    print 'Merging ligand with PDB Protein'

    clu = mcproPath + '/miscexec/clu '

    outputPDB = complexPDB[:-4] + '_cplx.pdb'

    if os.path.isfile(outputPDB): os.remove(outputPDB)

    cmd = clu + ' -t:s=' + resnumber + ' ' + complexPDB + ' -r ' + ligandZmat + ' -n ' + outputPDB

    print cmd

    os.system(cmd)

def generateScript(complexPDB, ligandLstToRemoveFromPDB, residueToBeTheCenterOfTheChopping, setCapOriginAtom,
                   setCutOffSize, HipLstOfResidues, HidLstOfResidues):

    lastPart = '''set variable origin ligand
set variable size 12
set targetq 0
fix charges
write pdb aaaa.chop.pdb
write pepz all aaaa.chop.all.in
write pepz variable aaaa.chop.var.in
write translation aaaa.chop.tt
EOF

'''

    complexPDBName = complexPDB[:-4]

    tmpScript = '$MCPROdir/miscexec/chop -u -i ' + complexPDBName + '_cplx.pdb << EOF\n'

    if len(ligandLstToRemoveFromPDB) != 0:
        for ele in ligandLstToRemoveFromPDB:
            tmpScript += 'delete ligand :' + ele + '\n'

    tmpScript += 'add center :' + residueToBeTheCenterOfTheChopping + '\n'
    tmpScript += 'set cap origin :' + residueToBeTheCenterOfTheChopping + '@' + setCapOriginAtom + '\n'

    tmpScript += 'set cut origin ligand \n'
    tmpScript += 'set cut size ' + setCutOffSize + '\n'
    #	tmpScript += 'delete cut :gol\n'

    tmpScript += 'fix chains\n'
#    tmpScript += 'delete cut :18c\n'
#    tmpScript += 'delete cut :22c\n'
    tmpScript += 'fix chains\n'
    tmpScript += 'cap all\n'
#    tmpScript += 'uncap :1a\n'
#    tmpScript += 'uncap :1b\n'
#    tmpScript += 'uncap :1c\n'


    if len(HipLstOfResidues) != 0:
        lst = ''
        for ele in HipLstOfResidues:
            lst += ele.strip() + ','
        tmpScript += 'set hip ' + ':' + lst[:-1] + '\n'

    if len(HidLstOfResidues) != 0:
        lst = ''
        for ele in HidLstOfResidues:
            lst += ele.strip() + ','
        tmpScript += 'set hid ' + ':' + lst[:-1] + '\n'

    tmpScript += lastPart.replace('aaaa', complexPDBName)

    return tmpScript

def prepareReducedChopped(complexPDB, ligandLstToRemoveFromPDB, residueToBeTheCenterOfTheChopping, setCapOriginAtom,
                          setCutOffSize, HipLstOfResidues, HidLstOfResidues):

    print 'Preparing Reduced Chopped System'

    chopScript = generateScript(complexPDB, ligandLstToRemoveFromPDB, residueToBeTheCenterOfTheChopping,
                                setCapOriginAtom, setCutOffSize, HipLstOfResidues, HidLstOfResidues)

    fileChopScript = open('chopScript.csh', 'w')
    fileChopScript.write(chopScript)
    fileChopScript.close()

    os.system('csh chopScript.csh')

def getNumberOfTheLastResidueOfTheChoppedSystem(choppedPDBName):

    resNumber = ''

    for line in open(choppedPDBName):
        if 'ATOM' in line[:6] or 'HETATM' in line[:6]:
            resNumber = line.split()[4]

    return resNumber

def prepareFinalZmatrixWithPEPz(complexPDB, titleOfYourSystem, ligandZmat, fixBackBoneSelection):

    complexPDBName = complexPDB[:-4]

    tmpfile = open(complexPDBName + '.all.in', 'w')

    for line in open(complexPDBName + '.chop.all.in'):
        newLine = line
        if '[ADD YOUR TITLE HERE]' in line: newLine = line.replace('[ADD YOUR TITLE HERE]', titleOfYourSystem)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: newLine = line.replace(
            '[WRITE NAME OF YOUR solute z-matrix FILE]', ligandZmat)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: newLine = line.replace(
            '[NAME OF THE z-matrix TO BE WRITTEN]', complexPDBName + 'all.z')

        tmpfile.write(newLine)

    tmpfile.close()

    tmpfile = open(complexPDBName + '.var.in', 'w')

    for line in open(complexPDBName + '.chop.var.in'):
        newLine = line
        if '[ADD YOUR TITLE HERE]' in line: newLine = line.replace('[ADD YOUR TITLE HERE]', titleOfYourSystem)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: newLine = line.replace(
            '[WRITE NAME OF YOUR solute z-matrix FILE]', ligandZmat)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: newLine = line.replace(
            '[NAME OF THE z-matrix TO BE WRITTEN]', complexPDBName + 'var.z')

        tmpfile.write(newLine)

    tmpfile.close()

    # create conrot

    lastResidue = getNumberOfTheLastResidueOfTheChoppedSystem(complexPDBName + '.chop.pdb')

    tmpfile = open(complexPDBName + '.var_conrot.in', 'w')

    for line in open(complexPDBName + '.var.in'):
        newLine = line
        if 'parameter type ALL *' in line: newLine = line + '$ set conrot\n$ set override domain 1-' + lastResidue + '\n'
        if '$ set fixed backbone' in line:
            newLine = ''
            for ele in fixBackBoneSelection:
                newLine += '$ set fixed backbone ' + ele + '\n'
        if '$ write zmatrix ' + complexPDBName + 'var.z' in line: newLine = line.replace(complexPDBName + 'var.z',
                                                                                         complexPDBName + 'varcon.z')
        tmpfile.write(newLine)

    tmpfile.close()

def fixZmatrix(matrixZ, ligandResidueName):

    print 'Fixing Z matrix TERZ .... ', matrixZ

    tmpfile = open('tmp.txt', 'w')

    be2allFile = [ele for ele in open(matrixZ)]

    for iter in range(len(be2allFile)):
        newLine = be2allFile[iter]
        if 'TERZ' in newLine:
            if ligandResidueName in be2allFile[iter + 1] or 'CAP' in be2allFile[iter + 1]:
                newLine = be2allFile[iter]
            else:
                continue

        tmpfile.write(newLine)

    tmpfile.close()

    os.system('cp tmp.txt ' + matrixZ)

def createZmatrices(complexPDB, ligandResidueName):

    complexPDBName = complexPDB[:-4]

    # delete the pepz zmat output files if exist
    try:
        os.remove(complexPDBName+'all.z')
        os.remove(complexPDBName+'var.z')
        os.remove(complexPDBName+'varcon.z')
    except OSError:
        pass

    os.system('xPEPZ ' + complexPDBName + '.all')
    os.system('xPEPZ ' + complexPDBName + '.var')
    os.system('xPEPZ ' + complexPDBName + '.var_conrot')

    fixZmatrix(complexPDBName + 'all.z', ligandResidueName)
    fixZmatrix(complexPDBName + 'var.z', ligandResidueName)
    fixZmatrix(complexPDBName + 'varcon.z', ligandResidueName)

def relaxProteinLigand(complexPDB, mcproPath):

    complexPDBName = complexPDB[:-4]
    try:
        os.mkdir('CG9')
    except:
        os.system('rm -r CG9')
        os.mkdir('CG9')

    os.system(
        'cp ' + complexPDBName + 'all.z CG9;cd CG9;' + mcproPath + '/scripts/xCGDD9 ' + complexPDBName + 'all')

def getOptimizedZmatrix():

    # By default the file is called optsum

    ZmatrixOptimized = []

    for line in open('CG9/optsum'):

        if 'Geometry Variations follow' in line: break

        ZmatrixOptimized.append(line)

    return ZmatrixOptimized

def replaceOptimizedCoordinatesInZmatrixFile(varFileName, capFileName, zmatrixOptimized):

    fileout = open(capFileName, 'w')

    # get bottom part of the var file

    read = False
    bottomDataInVarFile = []

    for line in open(varFileName):
        if 'Geometry Variations follow' in line: read = True
        if read: bottomDataInVarFile.append(line)

    for ele in zmatrixOptimized: fileout.write(ele)
    for ele in bottomDataInVarFile: fileout.write(ele)

def generateFinalStructuresWithCAP(complexPDB):

    print 'Generating final structures'

    complexPDBName = complexPDB[:-4]

    zmatrixOptimized = getOptimizedZmatrix()

    replaceOptimizedCoordinatesInZmatrixFile(complexPDBName + 'var.z', complexPDBName + 'cap.z', zmatrixOptimized)
    replaceOptimizedCoordinatesInZmatrixFile(complexPDBName + 'varcon.z', complexPDBName + 'capcon.z',
                                             zmatrixOptimized)

    try:
        os.mkdir('finalZmatrices')
    except:
        os.system('rm -r finalZmatrices')
        os.mkdir('finalZmatrices')

    os.system('cp CG9/optsum finalZmatrices/' + complexPDBName + 'all.z')
    os.system('cp ' + complexPDBName + 'cap.z ' + complexPDBName + 'capcon.z finalZmatrices')



        #  ********************************************************************************
if __name__ == "__main__":

	mcproPath = os.environ.get('MCPROdir')
	BOSSPath = os.environ.get('BOSSdir')

	BOSSscriptsPath = os.path.join(BOSSPath,'scripts')

	#  *********************** CODE STARTS *********************************************

checkParameters()

input = raw_input("Is optimized zmat available? y/n:"+'\n')
# input ='y'
if input=='y':
    print "Taking opted zmat from" ,ligandZmatOrig
    print "skipping GenLigPDB, PreptLigWChimera and OptLig"
    _,resnumber, resnumberLigandOriginal = generateLigandPDB(complexPDB,ligandResidueName)
    ligandZmat = fixDummyAtomNamesDirect(ligandZmatOrig)


else:
    print "Constructing zmat from the complex file"
    ligandZmat, resnumber, resnumberLigandOriginal = prepareLigandZmatrix(complexPDB, ligandResidueName, mcproPath,
                                                                  BOSSscriptsPath)

fixedComplexPDB = fixPDBprotein(complexPDB, ligandResidueName)

mergeLigandWithPDBProtein(mcproPath,fixedComplexPDB,ligandZmat,resnumberLigandOriginal)

prepareReducedChopped(fixedComplexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,HipLstOfResidues,HidLstOfResidues)

input = raw_input("\n Any additional modification on residue before executing pepz? y/n:"+'\n')
#
if input=='y':

    resnumberlist,reslinelist=locateRes(fixedComplexPDB, ResidueChangeSelection,ResidueFinalSelection)

    ReplaceRes(fixedComplexPDB,resnumberlist,reslinelist,ResidueFinalSelection)

FixShortChain(fixedComplexPDB,fixShortChainCri)

prepareFinalZmatrixWithPEPz(fixedComplexPDB,titleOfYourSystem,ligandZmat,fixBackBoneSelection)

createZmatrices(fixedComplexPDB, ligandResidueName)

relaxProteinLigand(fixedComplexPDB,mcproPath)

generateFinalStructuresWithCAP(fixedComplexPDB)
