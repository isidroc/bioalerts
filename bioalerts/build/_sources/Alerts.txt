Alerts: Derivation of structural alerts
*****************************

This module is composed of the following three classes:

* GetSubstructuresFromReferenceDataset
* CalculatePvaluesCategorical
* CalculatePvaluesContinuous

GetSubstructuresFromReferenceDataset
====================================
.. py:class:: GetSubstructuresFromReferenceDataset(name_field=None)

    Similar to GetDataSetInfo from the module LoadMolecules.
    
    This class can be used to load molecules in smiles, sdf and mol2 format.
    
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


CalculatePvaluesCategorical
===========================
.. py:class:: CalculatePvaluesCategorical(max_radius)

    :var max_radius: maximum substructure radius considered.
    :vartype max_radius: int
    :var output: pandas dataframe containing the following columns: (i) 'Compound ID', (ii) 'p_value', (iii) 'Compounds with substr.', (iv) 'activity label', (v) 'Substructure', (vi) 'Substructure in Molecule', (vii) 'Comp. with substr. active', and (viii) 'Comp. with substr. inactive'.
     :vartype output: pandas.core.frame.DataFrame


    .. py:method:: calculate_p_values(mols,substructure_dictionary,bioactivities,mols_ids,threshold_nb_substructures,threshold_pvalue,threshold_frequency,active_label=1,inactive_label=0)
    .. py:method:: HTMLOutputWriter(self,output_filename)
    .. py:method:: XlSXOutputWriter(self,frame, output_filename, molCol=['Substructure',"Substructure in Molecule"], size=(300,300))



CalculatePvaluesContinuous
==========================
.. py:class:: CalculatePvaluesContinuous(radii_ext)

    :var radii_ext: substructure radii to be considered
    :vartype radii_ext: list
    :var output: pandas dataframe containing the following columns: (i) 'Compound ID', (ii) 'Number compounds', (iii) 'statistic', (iv) 'p_value', (v) 'Diff. distribution means (w - wo)', (vi) 'Compounds with substr.', (vii) 'Substructure', and (viii) 'Substructure in Molecule'.
    :vartype output: pandas.core.frame.DataFrame


    .. py:method:: calculate_p_values(mols,substructure_dictionary,bioactivities,mols_ids,threshold_nb_substructures,threshold_pvalue,threshold_ratio_w_wo,Bonferroni=True)
    .. py:method:: HTMLOutputWriter(self,output_filename)
    .. py:method:: XlSXOutputWriter(self,frame, output_filename, molCol=['Substructure',"Substructure in Molecule"], size=(300,300))





