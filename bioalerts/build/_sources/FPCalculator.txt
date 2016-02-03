FPCalculator: Fingerprint Calculator 
*************************************

This module can be used to extract the substructures, of a user-defined radius,
from a list of molecules in rdkit Mol class format.

The module is composed of the following class:

* CalculateFPs

CalculateFPs
============

.. py:class:: CalculateFPs(radii,mols,reference_substructure_keys={})


    *CalculateFPs* serve to calculate the following types of Morgan fingerpints:

    1. Hashed fingerprints in binary format
    2. Hashed fingerprints in count format
    3. Unhashed fingerprints in binary format
    4. Unhashed fingerprints in count format

    Morgan fingerprints encode chemical structures by considering atom neighbourhoods. 
    Each substructure in the molecule set, with a maximal user-defined bond radius, is assigned an unambiguous integer identifier. 
    These identifiers are mapped either into an unhashed or hashed array. 
    For the hashed array, the position in the array where the substructures will be mapped is given by the modulo of the division 
    of the substructure identifier by the fingerprint size. 
    In the case of unhashed (keyed) fingerprints, each bit in the fingerprint is associated to only one substructure, 
    producing a length of the unhashed fingerprints equal to the number of distinct substructures present in the dataset.
    Both hashed and unhashed fingerprints can be stored in binary and count format. 
    In count format, each bit in the fingerprint accounts for the number of times each substructure is present in a given compound,
    whereas in binary format each bit encodes whether a substructure is present in a compound (1), irrespective of the number of occurrences, or not (0).


    :var radii: radii of the substructures that will be used to generate the fingerprints for the molecules specified in the argument *mols*.
    :var max_radius: maximum substructure radius considered.
    :var mols: input molecules for which fingerprints will be calculated.
    :var reference_substructure_keys: dictionary of substructure identifiers that will be used to calculate the unhashed version of the Morgan fingerprints. This dictionary is calculated with the method *LoadMolecules.GetDataSetInfo.extract_substructure_information().* 
    :var substructure_dictionary: dictionary containing substructure for the molecules specified in the argument *mols*, and for which fingerprints will be calculated. 
    :var mols_reference_for_unhashed: reference set of molecules whose substructures are to be considered when computing unhashed fingerprints. This molecule set can be the same set of molecules for which the user wants to compute the fingerprints, *i.e.* *mols*, or a different molecule set. For instance, if a bioactivity model trained on a given data set using unhashed fingerprints is to be applied on an external data set, the fingerprints for the new molecules should be the same as those used to train the model. Thus, the refernce set of molecules (*mols_reference_for_unhashed*) would correspond in this case to the molecules used to train the aforesaid model. 
    
    The value for the following attributes is set when running the fingerprint calculation methods explained below.

    :var columns_unhashed: fingerprint identifiers corresponding to the columns in the unhashed fingerprints.
    :var substructure_ids: identifiers for the substructures.
    :var fps_hashed_binary_quick: hashed fingerprints in binary format
    :var fps_hashed_binary:  hashed fingerprints in binary format.
    :var fps_hashed_counts: unhashed fingerprints in count format.
    :var fps_unhashed_binary: unhashed fingerprints in binary format.
    :var fps_unhashed_counts: unhashed fingerprints in count format.
    :var substructures_smiles: dictionary containing the smiles for the substructures.
    :vartype columns_unhashed: numpy.ndarray 
    :vartype substructure_ids: numpy.ndarray
    :vartype fps_hashed_binary_quick: numpy.ndarray
    :vartype fps_hashed_binary: numpy.ndarray
    :vartype fps_hashed_counts: numpy.ndarray
    :vartype fps_unhashed_binary: numpy.ndarray
    :vartype fps_unhashed_counts: numpy.ndarray
    :vartype substructures_smiles: dict


    .. py:method:: calculate_hashed_fps_binary_quick(nBits)
        
        Fast class method to compute hashed fingerprints in binary format. The computed hashed fingerprints are stored in *CalculateFPs.fps_hashed_binary_quick*.
        
        :param nBits: fingerprint size
        :type nBits: int
    
    .. py:method:: calculate_hashed_fps_counts(nBits)
    
        Class method to compute hashed fingerprints in binary and count format. The computed fingerprints are stored in *CalculateFPs.fps_hashed_binary* and *CalculateFPs.fps_hashed_counts*, respectively.
        
        :param nBits: fingerprint size
        :type nBits: int
    
    .. py:method:: calculate_unhashed_fps(draw_substructures=False,image_directory='./images_substructures')
    
        Class method to compute unhashed fingerprints in binary and count format. If *draw_substructures* is set to True, depictions of each substructure in .pdf format, in the context of a molecules from *mols* where the substructure is present, will be saved to the directory specified in the argument *image_directory*.
        
        :param draw_substructures: if set to True, depictions for all substructures (with a bond radius allowed by the user through the argument radii) present in the training set of molecules will be generated.
        :type draw_substructures: bool
        :param image_directory: directory to save the substructure depictions
        :type image_directory: str



