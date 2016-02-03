from numpy import savetxt as _savetxt
from rdkit.Chem import SmilesMolSupplier as _SmilesMolSupplier, MolFromMol2Block as _MolFromMol2Block, SDMolSupplier as _SDMolSupplier
from os.path import splitext as _splitext, exists as _exists 
from operator import add as _add
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect as  _GetMorganFingerprintAsBitVect, GetMorganFingerprint as _GetMorganFingerprint



def RetrieveMol2Block(fileLikeObject, delimiter="@<TRIPOS>MOLECULE"):
    """Generator to retrieve one mol2 block at a time when parsing a mol file
    """
    mol2 = []
    for line in fileLikeObject:
        if line.startswith(delimiter) and mol2:
            yield "".join(mol2)
            mol2 = []
        mol2.append(line)
    if mol2:
        yield "".join(mol2)


class LoadMolecules:    
    """Load molecules from (i) smiles, (ii) sdf, and (iii) mol2 files.
    """
    def __init__(self,input_file,verbose=True,delimiter="\t",name_field="_Name"):
        self.input_file = input_file
        self.verbose = True
        self.delimiter = delimiter
        file_name, file_extension = _splitext(input_file)
        self.file_name = file_name
        self.file_extension = file_extension
        if(file_extension not in ['.smi','.smiles','.sdf','.mol2']): 
            raise ValueError("Incorrect file extension")
        self.mols = []
        self.molserr = []
        self.nb_mols = None
        self.mols_ids = []
        self.name_field = name_field
        
    def ReadMolecules(self,titleLine=False,smilesColumn=0,nameColumn=1): #titleLine for smiles

        if self.file_extension in ['.smi','.smiles']:
            if self.verbose:
                print "Format of the structures file = SMILES"
            suppl = _SmilesMolSupplier(self.input_file,smilesColumn=smilesColumn,
                                           nameColumn=nameColumn,
                                           delimiter=self.delimiter,titleLine=titleLine)

            for i,m in enumerate(suppl):
                if m is not None:
                    self.mols.append(m)
                    mol_id = i if self.name_field == None else m.GetProp(self.name_field)
                    self.mols_ids.append(mol_id)
                else:
                    self.molserr.append(i)
            nb_mols=len(self.mols)
        elif self.file_extension == '.mol2':
            print "Format of the structures file = Mol2"
            molss=[]
            with open(self.input_file) as fi:
                for mol2 in RetrieveMol2Block(fi):
                    rdkMolecule = _MolFromMol2Block(mol2)
                    molss.append(rdkMolecule)
            for i,m in enumerate(molss):
                if m is not None:
                    self.mols.append(m)
                    mol_id = i if self.name_field == None else m.GetProp(self.name_field)
                    self.mols_ids.append(mol_id)
                else:
                    self.molserr.append(i)
                    self.mols.append(m)
            self.nb_mols=len(self.mols)
        else:
            if self.verbose:
                print "Format of the structures file = SDF"
            suppl = _SDMolSupplier(self.input_file)
            for i,m in enumerate(suppl):
                if m is not None:
                    self.mols.append(m)
                    mol_id = i if self.name_field == None else m.GetProp(self.name_field)
                    self.mols_ids.append(mol_id)
                else:
                    self.molserr.append(i)
            self.nbMols=len(self.mols)
        
        if self.verbose:
            if len(self.molserr) !=0:
                print "%d molecules (starting at zero) could not be processed.\n"%(len(self.molserr))
                err_file="incorrect_molecules.csv"
                print "This information has been saved in the following file: %s\n"%(err_file)
                #for x in self.molserr: print x
                print "NOTE: the indexes of the molecules start at zero. Thus the first molecule is molecule 0."
                # Save the information about which molecules could not be processed correctly.
                _savetxt(err_file,self.molserr,fmt="%d")
            else:
                print "All molecules in the input file were processed correctly"

class GetDataSetInfo:
    '''
    Crate a dictionary: keys = substructure IDs, value = compound IDs.
        Thus, we know for a compound, which substructures it contains
    '''

    def __init__(self,name_field=None):
        self.name_field = name_field
        self.nb_substructures = None
        self.max_radius = None
        self.mols_ids = []
        self.substructure_dictionary = {}
    
    def _combine_dicts(self,a, b, op=_add):
        return dict(a.items() + b.items() + [(k, op(a[k], b[k])) for k in set(b) & set(a)])
        
    def extract_substructure_information(self,radii,mols):
        self.radii = radii
        global indexes_mols        
        for i,m in enumerate(mols):
               info={}
               fp = _GetMorganFingerprint(m,max(radii),bitInfo=info)
               mol_id = i if self.name_field == None else m.GetProp(self.name_field)
               self.mols_ids.append(mol_id)
               substructure_dictionary = {k:[mol_id] for k,v in info.iteritems() if v[0][1] in radii}
               self.substructure_dictionary = self._combine_dicts(substructure_dictionary,self.substructure_dictionary)
        self.nb_substructures = len(self.substructure_dictionary.keys())
        self.max_radius = max(radii)



