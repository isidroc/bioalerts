#!/usr/bin/python

# Isidro Cortes Ciriano.    6/8/2013
# Institut Pasteur
# isidrolauscher@gmail.com

# Import modules
import argparse
import numpy as np
import os,sys
# Arguments passed to the scripts
parser = argparse.ArgumentParser(prog='PROG',description='Get Morgan Fingerprints for compounds codified in either SMILES or SDF format using RDkit. Isidro Cortes Ciriano. August/September 2013')
parser.add_argument('--bits', required='TRUE',type=int, help="Size of the hashed Morgan Fingerprints (binary and with counts)")
parser.add_argument('--rad', required='TRUE', type=int, help="Maximum radius of the substructures. Deafault is two, equivalent to ECFP4 from PipelinePilot")
parser.add_argument('--f', required='TRUE', type=str, help="Format of the input file")
parser.add_argument('--mols', type=str,help="File containing the molecules {.smi|.smiles|.sdf|.mol2}. If the format is smiles, each line should contain the smiles and the name separated by a comma (in this order)")
parser.add_argument('--image', action='store_true', help="Write --image if you want the images of the substructures")
parser.add_argument('--unhashed', action='store_true', help="Write --unhashed if you want the unhashed fingerprints")
parser.add_argument('--v', action='store_true', help="Verbose")
parser.add_argument('--extF',type=str, help="Type -extF followed by the format {.smi|.smiles|.sdf} of the external file for which you want to calculate HASHED circular fingerprints")
parser.add_argument('--molsEXT',type=str,help="External file")
parser.add_argument('--unhashedEXT', action='store_true', help="Write --unhashedEXT if you want the unhashed fingerprints for the external file. The substructures of the molecules in the external file will be compared to the pool of substructures contained in the molecules of the main file")
parser.add_argument('--RDkitPath', required='TRUE', type=str, help="Path to the directory where the RDkit files are")
parser.add_argument('--output', required='TRUE', type=str, help="Name of the output files")
args = vars(parser.parse_args())

image=args['image']
unhashed=args['unhashed']
verbose=args['v']
formatFile=args['f']
fileMols=str(args['mols'])
nbBits=int(args['bits'])
fp_diam=int(args['rad'])
# External file.
formatFileEXT=args['extF']
fileMolsEXT=str(args['molsEXT'])
unhashedEXT=args['unhashedEXT']
RDkitPath=args['f']
outname=args['output']
sys.path.append(RDkitPath)


#if (formatFileEXT and not fileMolsEXT) or (fileMolsEXT and not formatFileEXT):
#    sys.exit("If molsEXT is defined, the argument extF also needs to be defined and vice versa.\nThe calculation has stopped here.")


if verbose:
    if image:
        print "\nCalculation of Morgan Fingerprints with diameter %d hashed into a fingerprint size equal to %d.\nMolecules file: %s.\nImages for the chemical substructures will be created.\n" %(args['rad'],args['bits'],args['mols'])
    else :
        print "\nCalculation of Morgan Fingerprints with diameter %d hashed into a fingerprint size equal to %d.\nMolecules file: %s.\n NO Images for the chemical substructures will be created.\n" %(args['rad'],args['bits'],args['mols'])

#####################################
# Import Modules
#####################################
import gzip
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import rdkit.rdBase
from rdkit.Chem.MACCSkeys import GenMACCSKeys
from rdkit.Chem import AllChem
from rdkit import DataStructs 
from rdkit.DataStructs import BitVectToText
from rdkit.Chem import Draw

#####################################
# Define Functions
#####################################
# To search within sublists:
def insublist(item, list):
    for l in list:
        if np.array_equal(item,l):
            return True
    return False

# Define a function to create matrices of empty strings
def nans(shape, dtype=str):
    a = np.empty(shape, dtype)
    a[:]=""
    return a

#####################################
# Open File
#####################################
# Open the files where the fingerprints will be kept:
fp_hash_b=outname+"_hashed_binary.csv"
if os.path.exists(fp_hash_b):
    os.remove(fp_hash_b)
f_fp_bin=open(fp_hash_b,'w')

fp_hash_c=outname+"_hashed_counts.csv"
if os.path.exists(fp_hash_c):
    os.remove(fp_hash_c)
f_fp_counts=open(fp_hash_c,'w')

#####################################
# Read Molecules
#####################################
### Read Mol2 files
def RetrieveMol2Block(fileLikeObject, delimiter="@<TRIPOS>MOLECULE"):
    import rdkit.Chem
    """generator which retrieves one mol2 block at a time
    """
    mol2 = []
    for line in fileLikeObject:
        if line.startswith(delimiter) and mol2:
            yield "".join(mol2)
            mol2 = []
        mol2.append(line)
    if mol2:
        yield "".join(mol2)



if formatFile == 'smi' or formatFile == 'smiles':
    if verbose:
        print "Format of the main file = SMILES"
    suppl = Chem.SmilesMolSupplier(fileMols,smilesColumn=0,nameColumn=1,delimiter=',',titleLine=False)
    mols=[]
    molserr=[]
    for i,m in enumerate(suppl):
        if m is not None:
            mols.append(m)
        else:
            molserr.append(i)
    nbMols=len(mols)
elif formatFile == 'mol2':
    print "molecules in mol2 format\n"
    molss=[]
    with open(fileMols) as fi:
        for mol2 in RetrieveMol2Block(fi):
            rdkMolecule = rdkit.Chem.MolFromMol2Block(mol2)
            molss.append(rdkMolecule)
    molserr=[]
    mols=[]
    for i,m in enumerate(molss):
        if m is not None:
            mols.append(m)
        else:
            molserr.append(i)
            mols.append(m)  
    nbMols=len(mols)
    print nbMols
else:
    if verbose:
        print "Format of the main file = SDF"
    suppl = Chem.SDMolSupplier(fileMols)
    mols=[]
    molserr=[]
    for i,m in enumerate(suppl):
        if m is not None:
            mols.append(m)
        else:
            molserr.append(i)
    nbMols=len(mols)

if verbose: 
    if len(molserr) !=0:
        print "The following %d molecules (starting at zero) could not be processed:\n"%(len(molserr))
        for x in molserr: print x
        print "NOTE: the indexes of the molecules start at zero. Thus the first molecule is molecule 0."
        errfile="incorrect_molecules_"+outname+".csv"
        print "This information has been saved in the following file: %s\n"%(errfile)
        # Save the information about which molecules could not be processed correctly.
        np.savetxt(errfile,molserr,fmt="%d")
        del errfile
    else:
        print "All molecules in the input file were processed correctly"

###########################
# External File
##########################
if formatFileEXT:
    molserrEXT=[]
    molsEXT=[]
    if formatFileEXT == 'smi' or formatFileEXT == 'smiles':
        if verbose:
            print "Format of the external file = SMILES"
        supplEXT = Chem.SmilesMolSupplier(fileMolsEXT,smilesColumn=0,nameColumn=1,delimiter=',',titleLine=False)
        for i,m in enumerate(supplEXT):
            if m is not None:
                molsEXT.append(m)
            else:
                molserrEXT.append(i)
        nbMolsEXT=len(molsEXT)
    elif formatFile == 'mol2':
        molssEXT=[]
        with open(fileMolsEXT) as fi:
            for mol2 in RetrieveMol2Block(fi):
                rdkMolecule = rdkit.Chem.MolFromMol2
                molssEXT.append(rdkMolecule)
        for i,m in enumerate(molssEXT):
            if m is not None:
                molsEXT.append(m)
            else:
                molserrEXT.append(i)
                molsEXT.append(m)  
        nbMolsEXT=len(molsEXT)
    else:
        if verbose:
            print "Format of the external file = SDF"
        supplEXT = Chem.SDMolSupplier(fileMolsEXT)
        for i,m in enumerate(supplEXT):
            if m is not None:
                molsEXT.append(m)
            else:
                molserrEXT.append(i)
        nbMolsEXT=len(molsEXT)

if verbose and formatFileEXT: 
    if len(molserrEXT) !=0:
        print "The following %d molecules (starting at zero) from the EXTERNAL file could not be processed:\n"%(len(molserr))
        for x in molserrEXT: print x
        print "NOTE: the indexes of the molecules start at zero. Thus the first molecule is molecule 0."
        errfileEXT="incorrect_molecules_EXT_"+outname+".csv"
        print "This information has been saved in the following file: %s\n"%(errfileEXT)
        # Save the information about which molecules could not be processed correctly.
        np.savetxt(errfileEXT,molserrEXT,fmt="%d")
        del errfileEXT
    else:
        print "All molecules in the EXTERNAL file were processed correctly"

    if verbose:
        print 'Your molecules file has %d CORRECT molecules\n' % (len(mols))
        if formatFileEXT:
            print 'Your external file contains %d CORRECT molecules\n' % (len(molsEXT))

#declare the vector of zeros to know which positions have appeared
position_track=[0]*nbBits

# Define variable to keep the smiles and bit numbers
arr = [[],[]]
for i in range(nbBits):
    arr[0].append(i)

for i in range(nbBits):
        arr[1].append([''])

# Define variable to keep the smiles and bit numbers
arr2 = [[],[]]
for i in range(nbBits):
    arr2[0].append(i)

for i in range(nbBits):
    arr2[1].append([''])


# Define the list of lists containing for each compound, the features that it contains.
fps_by_comp=[[]]
for i in range(nbMols):
    fps_by_comp[0].append([''])

# Define the list that will contain all the submolecules
subm_all=[]
# Define the list where the erroneous molecules will be saved.
err_mols=[]

# Define a progess bar
from progressbar import *
widgets = ['Progression: ', Percentage(), ' ', Bar(marker='.',left='[',right=']'),' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

pbar = ProgressBar(widgets=widgets, maxval=nbMols)
pbar.start()

smiles_subs_kept =[]
Atoms_subs=[]
nbFeatTot=0

#Loop over the molecules
for molecule_nb,m in enumerate(mols):
    info={}; info2={}
    if m is None:
            print "Erroneous input at molecule: %d" %(molecule_nb)
            err_mols.append(molecule_nb)
    else:
        if image:
            image_name="%s_Molecule_%d.pdf"%(outname,molecule_nb+1)
            tmp=AllChem.Compute2DCoords(m)
            Draw.MolToFile(m,image_name,size=(300,300),wedgeBonds=True,kekulize=True)
#        if verbose:
#            print "Molecule %d\n" % (molecule_nb)
        fp = AllChem.GetMorganFingerprintAsBitVect(m,fp_diam,nbBits,bitInfo=info) 
        AllBits=np.asarray([info.items()[i][0] for i in range(0,len(info.items()))])
        #diameter 2 is equal to the length 4 of ECFP-4
        fp_bits=BitVectToText(fp)
        fp_counts=list(fp_bits)
        AtomRadBits=np.asarray([info.items()[p][1][:] for p in range(0,len(info.items()))])

        fp2 = AllChem.GetMorganFingerprint(m,fp_diam,bitInfo=info2)
        ids_now=np.asarray([info2.items()[i][0] for i in range(0,len(info2.items()))])
        ids_nowMOD=ids_now%nbBits
        AtomRadNow=np.asarray([info2.items()[p][1][:] for p in range(0,len(info2.items()))])

        for i in range(0,len(info2.items())):
            radius=info2.items()[i][1][0][1]
            atom=info2.items()[i][1][0][0]

            for k in range(0,len(AllBits)):
                if insublist(AtomRadNow[i][0],AtomRadBits[k]): 
                    bit = AllBits[k]
                    break

            counts=len(info2.items()[i][1])
            if counts > 1: 
                fp_counts[bit]=str(counts)

            if position_track[bit] == 0: 
                position_track[bit]=1
                env=Chem.FindAtomEnvironmentOfRadiusN(m,radius,atom)
                amap={}
                submol=Chem.PathToSubmol(m,env,atomMap=amap)
                if radius ==0: ##if len(amap)==0: # This means that the radius is zero, so the feature is a single atom
                    pass
                #    arr[1][bit].append(ids_now[i]) #submol
                #    # Draw the feature
                #    if image and ids_now[i] not in subm_all:
                #        image_name="%s_Feature_%d.pdf"%(outname,nbFeatTot)
                #        amap={}; amap[atom] = atom    
                #        Draw.MolToFile(m,image_name,size=(300,300),highlightAtoms=amap.keys())
                #    smiles_subs_kept.append(Chem.MolToSmiles(submol))
                #    Atoms_subs.append(submol.GetNumAtoms())
                #    if ids_now[i]  not in subm_all:
                #        nbFeatTot+=1 
                #        subm_all.append(ids_now[i])
                #    # For each molecule keep the substructures
                #    fps_by_comp[0][molecule_nb].append(ids_now[i])
                #    arr2[1][bit].append(str(nbFeatTot))
                else:
                    arr[1][bit].append(ids_now[i] )
                    if image and ids_now[i]  not in subm_all:
                        image_name="%s_Feature_%d.pdf"%(outname,nbFeatTot)
                        Draw.MolToFile(m,image_name,size=(300,300),highlightAtoms=amap.keys())
                    smiles_subs_kept.append(Chem.MolToSmiles(submol))
                    Atoms_subs.append(submol.GetNumAtoms())
                    if ids_now[i] not in subm_all:
                        nbFeatTot+=1
                        subm_all.append(ids_now[i] )
                    fps_by_comp[0][molecule_nb].append(ids_now[i] ) 
                    arr2[1][bit].append(str(nbFeatTot)) 


            else: #The bit is already on!
                env=Chem.FindAtomEnvironmentOfRadiusN(m,radius,atom)
                amap={}
                submol=Chem.PathToSubmol(m,env,atomMap=amap)

                if radius == 0 and ids_now[i]  not in arr[1][bit]:  ####:if len(amap)==0 and ids_now[i]  not in arr[1][bit]:
                    pass
                #    arr[1][bit].append(ids_now[i] )
                #    if image and ids_now[i]  not in subm_all:
                #        image_name="%s_Feature_%d.pdf"%(outname,nbFeatTot)
                #        amap={}; amap[atom] = atom
                #        Draw.MolToFile(m,image_name,size=(300,300),highlightAtoms=amap.keys())
                #    smiles_subs_kept.append(Chem.MolToSmiles(submol))
                #    Atoms_subs.append(submol.GetNumAtoms())
                #    if ids_now[i]  not in subm_all:
                #        nbFeatTot+=1
                #        subm_all.append(ids_now[i])
                #    fps_by_comp[0][molecule_nb].append(ids_now[i])
                #    arr2[1][bit].append(str(nbFeatTot)) 
                # We keep the all the features for each compound anyway
                if radius == 123123 and ids_now[i]  in arr[1][bit]: ###if len(amap)==0 and ids_now[i]  in arr[1][bit]:
                    pass
                    #fps_by_comp[0][molecule_nb].append(ids_now[i] )

                if submol.GetNumAtoms() >1 and ids_now[i]  not in arr[1][bit]: #####len(amap)!=0 and ids_now[i]  not in arr[1][bit]:
                    arr[1][bit].append(ids_now[i] )
                    if image and ids_now[i]  not in subm_all:
                        image_name="%s_Feature_%d.pdf"%(outname,nbFeatTot)
                        Draw.MolToFile(m,image_name,size=(300,300),highlightAtoms=amap.keys())
                    smiles_subs_kept.append(Chem.MolToSmiles(submol))
                    Atoms_subs.append(submol.GetNumAtoms())
                    if ids_now[i] not in subm_all:
                        nbFeatTot+=1
                        subm_all.append(ids_now[i] )
                    fps_by_comp[0][molecule_nb].append(ids_now[i] ) 
                    arr2[1][bit].append(str(nbFeatTot)) 
                # We keep the all the features for each compound anyway
                if submol.GetNumAtoms() >1 and ids_now[i]  in arr[1][bit]: ####if len(amap)!=0 and ids_now[i]  in arr[1][bit]:
                    fps_by_comp[0][molecule_nb].append(ids_now[i] )

         # Print the features in the corresponding files
        count=1
        for item in fp_bits:
             if count != nbBits:
               f_fp_bin.write("%s," % (item))
             else:
               f_fp_bin.write("%s" % (item))
             count+=1
        f_fp_bin.write("\n")
           
        count=1
        for item in fp_counts:
           if count != nbBits:
             f_fp_counts.write("%s," % (item))
           else:
             f_fp_counts.write("%s" % (item))
           count+=1
        f_fp_counts.write("\n")

# Updating the progress bar.
    #if verbose:
    if nbMols % (1+molecule_nb) == 0:
        pbar.update(molecule_nb)
        print "\n"

f_fp_bin.close()
f_fp_counts.close()

fp_per_bit=outname+'_features_per_bit_hashed_fp.csv'
if os.path.exists(fp_per_bit):
    os.remove(fp_per_bit)


f=open(fp_per_bit,'w')
for i in arr2[0]:
     if arr2[1][i]==['']:
         f.write("%d\n" % (i))
     else:
      f.write("%d%s\n" % (i, ','.join(list(set(arr2[1][i])))))

f.close()


if verbose:
    print "Total number of features : %d" %(len(subm_all))

if unhashed:
    FPS=nans((nbMols,len(subm_all)))
    FPS_counts=nans((nbMols,len(subm_all)))
    if verbose:
        print "Writing UNhashed fingerprints to file.."
    for i in range(nbMols):
        for j in range(len(subm_all)):
            if subm_all[j] in fps_by_comp[0][i]:
                FPS[i][j]=1
                FPS_counts[i][j]=fps_by_comp[0][i].count(subm_all[j])
            else:
                FPS[i][j]=0
                FPS_counts[i][j]=0

    fpbinary=outname+'_unhashed_binary.csv'
    fpcounts=outname+'_unhashed_counts.csv'
    np.savetxt(fpbinary, FPS, fmt='%1s', delimiter=',', newline='\n')
    np.savetxt(fpcounts, FPS_counts, fmt='%1s', delimiter=',', newline='\n')



###############################
# Write the smiles for the substructures
###############################

filename = outname+"_smiles_substructures.smi"
f = open(filename,'w')
dat = 'Substructure_ID\tNum_Atoms\tSmiles\n'
f.write(dat)
for i,m in enumerate(smiles_subs_kept):
    dat = i+'\t'+str(Atoms_subs[i])+'\t'+m+'\n'
    f.write(dat)
f.close()
###############################
# External Dataset
###############################
if formatFileEXT:
#    if verbose:
#        print "\nProcessing the external file..\n"
###############################
# Open File
###############################
# Open the files where the fingerprints will be kept:
    binaryEXT=outname+"_hashed_binary_EXT.csv"
    countsEXT=outname+"_hashed_counts_EXT.csv"
    if os.path.exists(binaryEXT):
        os.remove(binaryEXT)
    f_fp_binEXT=open(binaryEXT,'w')

    if os.path.exists(countsEXT):
        os.remove(countsEXT)
    f_fp_countsEXT=open(countsEXT,'w')


#declare the vector of zeros to know which positions have appeared
    position_trackEXT=[0]*nbBits

# Define variable to keep the smiles and bit numbers
    arrEXT = [[],[[]]]
    for i in range(nbBits):
        arrEXT[0].append(i)

    for i in range(nbBits):
            arrEXT[1].append([''])



# Define the list of lists containing for each compound, the features that it contains.
    fps_by_compEXT=[[]]
    for i in range(nbMols):
        fps_by_compEXT[0].append([''])

#Loop over the molecules
    for molecule_nb,m in enumerate(molsEXT):
        infoFP={}; infoEXT={}
        if m is None:
                print "Erroneous input at molecule (external file): %d" %(molecule_nb)
        else:
            if image:
                image_name="Molecule_Ext_%d.pdf"%(molecule_nb+1)
                tmp=AllChem.Compute2DCoords(m)
                Draw.MolToFile(m,image_name,size=(300,300),wedgeBonds=True,kekulize=True)
            #if verbose:
            #    print "External molecule: %d\n" % (molecule_nb)
            fpEXT = AllChem.GetMorganFingerprintAsBitVect(m,fp_diam,nbBits,bitInfo=infoFP) 
            fp_bitsEXT=BitVectToText(fpEXT)
            fp_countsEXT=list(fp_bitsEXT)
            AllBitsEXT=np.asarray([infoFP.items()[i][0] for i in range(0,len(infoFP.items()))])
            AtomRadBitsEXT=np.asarray([infoFP.items()[p][1][:] for p in range(0,len(infoFP.items()))])

            fp2EXT = AllChem.GetMorganFingerprint(m,fp_diam,bitInfo=infoEXT)
            ids_nowEXT=np.asarray([infoEXT.items()[i][0] for i in range(0,len(infoEXT.items()))])
            ids_nowMODEXT=ids_nowEXT%nbBits
            AtomRadNowEXT=np.asarray([infoEXT.items()[p][1][:] for p in range(0,len(infoEXT.items()))])

            for i in range(0,len(infoEXT.items())):
                radius=infoEXT.items()[i][1][0][1]
                atom=infoEXT.items()[i][1][0][0]

            for k in range(0,len(AllBitsEXT)):
                if insublist(AtomRadNowEXT[i][0],AtomRadBitsEXT[k]): 
                    bit = AllBitsEXT[k]
                    break

            countsEXT=len(infoEXT.items()[i][1])
            if countsEXT > 1: 
                fp_countsEXT[bit]=str(countsEXT)
            if position_trackEXT[bit] == 0: 
                position_trackEXT[bit]=1
            arrEXT[1][bit].append(ids_nowEXT[i]) #submol
            fps_by_compEXT[0][molecule_nb].append(ids_nowEXT[i] )

        # Writing the EXTERNAL fingerprints
        count=1
        for item in fp_bitsEXT:
             if count != nbBits:
               f_fp_binEXT.write("%s," % (item))
             else:
               f_fp_binEXT.write("%s" % (item))
             count+=1
        f_fp_binEXT.write("\n")
           
        count=1
        for item in fp_countsEXT:
           if count != nbBits:
             f_fp_countsEXT.write("%s," % (item))
           else:
             f_fp_countsEXT.write("%s" % (item))
           count+=1
        f_fp_countsEXT.write("\n")

    f_fp_binEXT.close()
    f_fp_countsEXT.close()


    if unhashedEXT:
        FPS_EXT=nans((nbMols,len(subm_all)))
        FPS_countsEXT=nans((nbMols,len(subm_all)))
        if verbose:
            print "Writing UNhashed fingerprints for the external file.."
        for i in range(nbMols):
            for j in range(len(subm_all)):
                if subm_all[j] in fps_by_compEXT[0][i]:
                    FPS_EXT[i][j]=1
                    FPS_countsEXT[i][j]=fps_by_compEXT[0][i].count(subm_all[j])
                else:
                    FPS_EXT[i][j]=0
                    FPS_countsEXT[i][j]=0

        outEXTbinary=outname+"_unhashed_binary_EXT.csv"
        outEXTcounts=outname+"_unhashed_counts_EXT.csv"
#        np.save("unhashed_binary.npy",FPS_EXT)
#        np.save("unhashed_counts.npy",FPS_countsEXT)
        np.savetxt(outEXTbinary, FPS_EXT, fmt='%1s', delimiter=',', newline='\n')
        np.savetxt(outEXTcounts, FPS_countsEXT, fmt='%1s', delimiter=',', newline='\n')

if verbose:
    print "Calculation Finished. NO problems encountered"

