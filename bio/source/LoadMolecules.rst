LoadMolecules: Load Molecules 
******************************

This module serves to load molecules from smiles, sdf or mol2 files.

The module is composed of the classes:

* LoadMolecules
* GetDataSetInfo

LoadMolecules
=============


.. py:class:: LoadMolecules (input_file,verbose=True,delimiter="\t",name_field="_Name")

This class can be used to load molecules to python from files in smiles, sdf and mol2 format
using the rdkit library (www.rdkit.org).
Molecules are thus saved as RDkit Mol class objects.
It contains the following instance attributes:

:var input_file: input file containing the molecules.
:vartype input_file: str
:var verbose: print to stdout the information about the progress of the calculations if set to True (default).
:vartype input_file: bool
:var delimiter: delimiter between fields in the input molecule file. Used with smiles files.
:vartype input_file: str
:var name_field: field in sdf files containing that will be used to get the molecule names.
:vartype name_field: str

The value for the following attributes are set when calling the method LoadMolecules.ReadMolecules():

:var mols: molecules that were properly processed.
:vartype mols: list
:var molserr: indices for the input molecules (starting at zero) that were incorrectly processed. This indices are also written to the file "incorrect_molecules.csv".
:vartype molserr: list
:var nb_mols: number of molecules properly processed.
:vartype nb_mols: int
:var mols_ids: indices of the molcules properly processed with respect to their position in the input file.
:vartype mols_ids: list
:var name_field: name of the field corresponding to the molecule names in the input file.
:vartype name_field: str


.. py:method:: ReadMolecules(titleLine=False,smilesColumn=0,nameColumn=1)

Method to read the molecules contained in the input file (LoadMolecules.input_file).
The arguments of this method are only used when reading molecules from smiles files. 

:param titleLine: From the RDkit documentation: "If true, the first line is assumed to list the names of properties in order seperated by 'delimiter'". Here, 'delimiter' corresponds to LoadMolecules.delimiter. 
:type titleLine: bool
:param smilesColumn: column (starting at zero) containing the smiles in the input file.
:type smilesColumn: int
:param nameColumn: column (starting at zero) containing the molecule names in the input file.
:type nameColumn: int

GetDataSetInfo
==============

.. py:class:: GetDataSetInfo(name_field=None)

This class can be used to load molecules in smiles sdf and mol2 format.

:var name_field: name of the data field containing the name of the molecules. Although the default value is None (*bool*), the name of the field if set when instantiating the class would be a string (*str*).
:vartype name_field: str
:var nb_substructures: total number of substructures, with a radius comprised in the argument *radii* of *GetDataSetInfo.extract_substructure_information()*, from the molecules specified in the argument *mols* of the method *GetDataSetInfo.extract_substructure_information()* (see below).
:vartype nb_substructures: int
:var max_radius: maximum substructure radius considered. This corresponds to the maximum value of the argument *radii* of *GetDataSetInfo.extract_substructure_information()*.
:vartype max_radius: int
:var substructure_dictionary: dictionary containing the substructures, with a radius comprised in the argument *radii* of *GetDataSetInfo.extract_substructure_information()*, and the molecules where they appear. Keys correspond to molecules names, whereas values correspond to 
:vartype substructure_dictionary: dict 

*GetDataSetInfo* contains the following methods:

.. py:method:: extract_substructure_information(radii,mols)
 
 This method extracts the substructures from the molecules (argument *mols*) 
whose radius is comprised in the argument *radii*.
Each substructure in the molecule set (*mols*) is assigned an unambiguous integer identifier,
which are kept in *GetDataSetInfo.substructure_dictionary*. 
The information about the substructures is kept in the fields indicated above.

:param radii: substructure radii to be considered
:param mols: molecules from which the substructures are to be extracted.
:type radii: list
:type mols: list


