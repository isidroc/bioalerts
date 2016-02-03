from numpy import arange as _arange, asarray as  _asarray, array as _array, sort as _sort, ones as _ones, zeros as  _zeros
from rdkit.Chem import SmilesMolSupplier as _SmilesMolSupplier, MolFromMol2Block as _MolFromMol2Block, SDMolSupplier as _SDMolSupplier, MolToSmiles as _MolToSmiles, PathToSubmol as _PathToSubmol, FindAtomEnvironmentOfRadiusN as _FindAtomEnvironmentOfRadiusN
from operator import add as _add
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect as  _GetMorganFingerprintAsBitVect, GetMorganFingerprint as _GetMorganFingerprint
from os.path import splitext as _splitext, exists as _exists 
from os import makedirs as _makedirs
from rdkit.Chem.Draw import MolToFile as _MolToFile

class CalculateFPs:
    '''
    Calculate fps for a set of molecules
    Required input: substructure dictionary, radii to be considered
    '''
    def __init__(self,radii,mols,reference_substructure_keys={}):
        self.radii = radii
        self.max_radius = max(radii)
        if type(mols) != list: mols = [ext.mols[i] for i in _arange(0,len(mols))] 
        self.mols = mols
        self.reference_substructure_keys = reference_substructure_keys
        self.substructure_dictionary = {}
        self.mols_reference_for_unhashed = None
        self.columns_unhashed = None 
        self.substructure_ids = None
        # output
        self.fps_hashed_binary_quick = None
        self.fps_hashed_binary = None
        self.fps_hashed_counts = None
        self.fps_unhashed_binary = None
        self.fps_unhashed_counts = None
        self.substructures_smiles = {}
        
    def _combine_dicts(self,a, b, op=_add):
        return dict(a.items() + b.items() + [(k, op(a[k], b[k])) for k in set(b) & set(a)])
    
    def calculate_hashed_fps_binary_quick(self,nBits):
        # bit format
        self.fps_hashed_binary_quick = _asarray([_GetMorganFingerprintAsBitVect(x,radius=self.max_radius,nBits=nBits) for x in self.mols])
            
    def calculate_hashed_fps(self,nBits):
        # count format
        fps_hashed_binary = _zeros((len(self.mols),nBits), dtype=int)
        fps_hashed_counts = _zeros((len(self.mols),nBits), dtype=int)
        for mol_index,mol in enumerate(self.mols): 
            info={}
            fp = _GetMorganFingerprint(mol,radius=self.max_radius,bitInfo=info)
            for key,val in info.iteritems():
                if val[0][1] in self.radii: #check if the radius is in the selection
                    fps_hashed_binary[mol_index,key%nBits] = 1
                    fps_hashed_counts[mol_index,key%nBits] += len(val)
        self.fps_hashed_binary = fps_hashed_binary
        self.fps_hashed_counts = fps_hashed_counts
                    
    
    def calculate_unhashed_fps(self,draw_substructures=False,image_directory='./images_substructures'): 
        # get the dictionary for the substructures
        idxs = []
        substr_ids = []
        counts = []    
        for mol_index,mol in enumerate(self.mols):
            info={}
            fp = _GetMorganFingerprint(mol,radius=self.max_radius,bitInfo=info)
            substructure_dictionary = {k:[mol_index] for k,v in info.iteritems() if v[0][1] in self.radii}
            substr_ids.append(substructure_dictionary.keys())
            idxs.append([mol_index]*len(substructure_dictionary.keys()))
            counts.append([ len(info.values()[x]) for x in _arange(0,len(info)) if info.values()[x][0][1] in self.radii])
            
            # get the smiles for the substructures
            amap = {}
            substructures_smiles = {k:[_MolToSmiles(_PathToSubmol(mol,_FindAtomEnvironmentOfRadiusN(mol,v[0][1],v[0][0]),atomMap=amap))] for k,v in info.iteritems() if v[0][1] in self.radii}
            self.substructures_smiles.update(substructures_smiles)
            
            # generate the images for the substructures if required..
            if draw_substructures:
                if not _exists(image_directory):
                    _makedirs(image_directory)
                for k,v in info.iteritems():
                    if k not in self.substructure_dictionary.keys() and v[0][1] in self.radii:
                        image_name="%s/Molecule_%d_substr_%d.pdf"%(image_directory,mol_index,k)
                        env=_FindAtomEnvironmentOfRadiusN(mol,v[0][1],v[0][0])
                        amap={}
                        submol=_PathToSubmol(mol,env,atomMap=amap)
                        _MolToFile(mol,image_name,size=(300,300),wedgeBonds=True,kekulize=True,highlightAtoms=amap.keys())
            
            self.substructure_dictionary = self._combine_dicts(substructure_dictionary,self.substructure_dictionary)
      
            
        idxs = _array([val for sublist in idxs for val in sublist])
        counts = _array([val for sublist in counts for val in sublist])
        substr_ids_flattened = [val for sublist in substr_ids for val in sublist]
        substr_ids = _array(substr_ids_flattened)
        self.substructure_ids = substr_ids
        if len(self.reference_substructure_keys)==0:
            print "No input set of keys for the substructures. \nThus, the substructures present in the input molecules will be considered for the calculation of unhashed fingerprints."
            columns = _array(list(set(self.substructure_dictionary.keys())))
            columns = _sort(columns)
            self.columns_unhashed = columns
            dimensionality_unhashed = len(columns)
        else:
            columns = _array(list(set(self.reference_substructure_keys)))
            columns = _sort(columns)
            self.columns_unhashed = columns
            dimensionality_unhashed = len(columns)
        
        fps_unhashed_binary = _zeros((len(self.mols),dimensionality_unhashed), dtype=int)
        fps_unhashed_counts = _zeros((len(self.mols),dimensionality_unhashed), dtype=int)
        

        # removing the indices corresponding to the substructures in the test molecules not present in the references set of substructures..
        idxs = _array([idxs[x] for x in _arange(0,len(substr_ids)) if substr_ids[x] in self.columns_unhashed])    
        counts = _array([counts[x] for x in _arange(0,len(substr_ids)) if substr_ids[x] in self.columns_unhashed]) 
        substr_ids = _array([substr_ids[x] for x in _arange(0,len(substr_ids)) if substr_ids[x] in self.columns_unhashed])
        mapping = _array([(substr_ids[x]==columns).nonzero() for x in _arange(0,len(substr_ids))])
        mapping = mapping.flatten()
        if len(mapping) ==0:
            print "There is no intersection between the substructures \n(i)provided in the reference key set, and\n(ii) the substructures found in the input molecules."
            return

        fps_unhashed_binary[idxs,mapping] = _ones(len(counts))
        fps_unhashed_counts[idxs,mapping] = counts
        self.fps_unhashed_binary = fps_unhashed_binary
        self.fps_unhashed_counts = fps_unhashed_counts
