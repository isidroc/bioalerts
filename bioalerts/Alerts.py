import numpy as np
import os
import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import rdkit.rdBase
from rdkit.Chem.MACCSkeys import GenMACCSKeys
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.DataStructs import BitVectToText
import scipy as sc
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from scipy import stats
from decimal import Decimal
import operator


class CalculatePvaluesCategorical:
    '''
    Calculate the p.values for each substructure in the input molecule
    n mols in the dataset
    m mols with a given activity label
    n' compounds with the substructure
    m' compounds from n' with the activity label
    '''
    def __init__(self,max_radius):
        # needs to be the same as the one used to calcualte the dictionary of substructures
        self.max_radius = max_radius  
        columns = ['Compound ID',
                   'Activity label',
                   'Substructure',
                   'Substructure in Molecule',
                   'p_value',
                   'Compounds with substr.',
                   'Comp. with substr. active',
                   'Comp. with substr. inactive',]
        self.output = pd.DataFrame(columns=columns)
        self.Bonferroni = None
    def calculate_p_values(self,mols,substructure_dictionary,bioactivities,mols_ids,threshold_frequency,
                           threshold_nb_substructures = 5,
                           threshold_pvalue = 0.05,
                           active_label=1,
                           inactive_label=0,
                           Bonferroni = True):
        self.Bonferroni = Bonferroni
        
        # n
        nb_mols = float(len(set([item for sublist in substructure_dictionary.values() for item in sublist])))
        # m
        nb_active_mols = float(np.sum(bioactivities == active_label))
        # (m - n)
        nb_inactive_mols = float(np.sum(bioactivities == inactive_label))

        nb_substructures_processed = 0
        if type(mols) != list: mols = [ext.mols[i] for i in np.arange(0,len(mols))]  #[x for x in mols]

        subs_discarded = [] # substructure that have been identified in other molecules. 
        for m,mol in enumerate(mols): #np.arange(0,len(mols)):
            #mol=mols[m]
            root_atoms_discarded = [] # center (or root) atoms discarded..
            info={}
            fp = AllChem.GetMorganFingerprint(mol,self.max_radius,bitInfo=info)
            # sort info to make sure the substructures are read from the smallest to the biggest.
            # In case a substructure with low radius is removed, we make sure all containing it will not be considered either in the following steps)
            # get keys sorted
            ff= sorted(info.iteritems(), key=operator.itemgetter(1))
            substructure_ids = [ff[x][0] for x in range(0,len(info))] 
            substructures_sub_dict = substructure_dictionary.keys()  

            for substructure_id in substructure_ids: 
                atom_radius = info[substructure_id]
                nb_substructures_processed += 1
                # check is the substructure is in the database (i.e. training data)
                if substructure_id in substructures_sub_dict and substructure_id not in subs_discarded and atom_radius[0][0] not in root_atoms_discarded:
                    mols_with_current_substructure = substructure_dictionary[substructure_id]
                    nb_comp_with_substructure = float(len(mols_with_current_substructure)) 
                    active_comp = (bioactivities == active_label)
                    comp_with_substructure = np.in1d(np.asarray(mols_ids) , np.asarray(mols_with_current_substructure))
                    nb_comp_with_substructure_active = np.sum(active_comp * comp_with_substructure) #i.e. m_{S act}
                    inactive_comp = (bioactivities == inactive_label)
                    #comp_with_substructure = np.in1d(np.asarray(mols_ids) , np.asarray(mols_with_current_substructure))
                    nb_comp_with_substructure_inactive = np.sum(inactive_comp * comp_with_substructure)

                    ## ACTIVE 
                    #########
                    #filter threshold of compounds with the substructure
                    filter_a = nb_comp_with_substructure > threshold_nb_substructures 
                    if filter_a: 
                        # filter threshold
                        filter_b = (float(nb_comp_with_substructure_active) / float(np.sum(comp_with_substructure))) > threshold_frequency
                        if filter_b:
                            p_value = 0 
                            for count in np.arange(nb_comp_with_substructure_active,nb_comp_with_substructure):
                                numerator = Decimal(sc.math.factorial(float(nb_comp_with_substructure)))
                                denominatorA = Decimal(sc.math.factorial(float(count))) * Decimal(sc.math.factorial(float(nb_comp_with_substructure-count)))
                                denominatorB = (nb_active_mols/nb_mols)**float(count)
                                denominatorC = (1.0 - (nb_active_mols/nb_mols))**(nb_comp_with_substructure - count)
                                out = float(numerator/denominatorA) * denominatorB * denominatorC
                                p_value += out 
                            
                            if p_value < threshold_pvalue:
                                #self.p_values_dictionary[substructure_id] = p_value
                        
                                # Drawing
                                env = Chem.FindAtomEnvironmentOfRadiusN(mol,atom_radius[0][1],atom_radius[0][0])
                                amap = {}
                                submol=Chem.PathToSubmol(mol,env,atomMap=amap)
                                m1 = mol
                                m1.GetSubstructMatch(submol)
                                #mm = Draw.MolToImage( mol,wedgeBonds=True,kekulize=True,highlightAtoms=amap.keys(),colour='green')
                                self.output = self.output.append({'Compound ID' : mols_ids[m],
                                                                  'Compounds with substr.' : nb_comp_with_substructure,
                                                                  #'Compounds with substr. and activity' : nb_comp_with_substructure_active,
                                                                  'p_value' : p_value,
                                                                  'Activity label':active_label,
                                                                  'Substructure in Molecule': m1,
                                                                  'Substructure':submol,
                                                                  'Comp. with substr. active':nb_comp_with_substructure_active,
                                                                  'Comp. with substr. inactive':nb_comp_with_substructure_inactive
                                                                  #'Smiles': Chem.MolToSmiles(mol)
                                                                  },ignore_index=True)
                                root_atoms_discarded.append(atom_radius[0][0])
                                subs_discarded.append(substructure_id)

                        
                        ## INACTIVE 
                        #########
                        #filter threshold of compounds with the substructure
                        # filter threshold
                        filter_b = (float(nb_comp_with_substructure_inactive) / float(np.sum(comp_with_substructure))) > threshold_frequency
                        if filter_b: 
                            p_value = 0 
                            for count in np.arange(nb_comp_with_substructure_inactive,nb_comp_with_substructure):
                                numerator = Decimal(sc.math.factorial(float(nb_comp_with_substructure)))
                                denominatorA = Decimal(sc.math.factorial(float(count))) * Decimal(sc.math.factorial(float(nb_comp_with_substructure-count)))
                                denominatorB = (nb_inactive_mols/nb_mols)**float(count)
                                denominatorC = (1.0 - (nb_inactive_mols/nb_mols))**(nb_comp_with_substructure - count)
                                out = float(numerator/denominatorA) * denominatorB * denominatorC
                                p_value += out 
                            
                            if p_value < threshold_pvalue:
                                #self.p_values_dictionary[substructure_id] = p_value
                        
                                # Drawing
                                env = Chem.FindAtomEnvironmentOfRadiusN(mol,atom_radius[0][1],atom_radius[0][0])
                                amap = {}
                                submol=Chem.PathToSubmol(mol,env,atomMap=amap)
                                m1 = mol
                                m1.GetSubstructMatch(submol)
                                #mm = Draw.MolToImage(mol,wedgeBonds=True,kekulize=True,highlightAtoms=amap.keys(),colour='red')
                                self.output = self.output.append({'Compound ID' : mols_ids[m],
                                                                  'Compounds with substr.' : nb_comp_with_substructure,
                                                                  #'Compounds with substr. and activity' : nb_comp_with_substructure_active,
                                                                  'p_value' : p_value,
                                                                  'Activity label':inactive_label,
                                                                  'Substructure in Molecule': m1,
                                                                  'Substructure':submol,
                                                                  'Comp. with substr. active':nb_comp_with_substructure_active,
                                                                  'Comp. with substr. inactive':nb_comp_with_substructure_inactive
                                                                  #'Smiles': Chem.MolToSmiles(mol)
                                                                  },ignore_index=True)
                                root_atoms_discarded.append(atom_radius[0][0])
                                subs_discarded.append(substructure_id)
                    else:
                        subs_discarded.append(substructure_id) 
                        root_atoms_discarded.append(atom_radius[0][0])
                                
        if self.Bonferroni == True:
           self.output['p_value'] = self.output['p_value'] * self.output.shape[0]
           self.output = self.output[self.output.p_value < 0.05]
        print 'Number of substructures processed: ', nb_substructures_processed 
        print 'Significant substructures: ', self.output.shape[0], 'substructures'

    def XlSXOutputWriter(self,output, output_filename, molCol=['Substructure',"Substructure in Molecule"], size=(300,300)):
        """
        http://nbviewer.ipython.org/github/Team-SKI/snippets/blob/master/IPython/rdkit_hackaton/XLSX%20export.ipynb
        Saves pandas DataFrame as a xlsx file with embedded images.
        It maps numpy data types to excel cell types:
        int, float -> number
        datetime -> datetime
        object -> string (limited to 32k character - xlsx limitations)
        
        Due to xlsxwriter limitations (other python xlsx writers have the same problem) 
        temporary image files have to be written to the hard drive.
        """
        
        import xlsxwriter 
           
        cols = list(output.columns)
        for i in molCol:
         cols.remove(i)

        dataTypes = dict(output.dtypes)
        
        workbook = xlsxwriter.Workbook(output_filename) # New workbook
        worksheet = workbook.add_worksheet() # New work sheet
        worksheet.set_column('A:A', size[0]/6.) # column width
        
        # Write first row with column names
        columns_nb = 0
        for x in molCol+cols:
            worksheet.write_string(0, columns_nb, x)
            columns_nb += 1
        
        row_nb = 1
        tmpfiles = []
        for index, row in output.iterrows():
            for j in np.arange(0,len(molCol)):
                imfile = "xlsx_tmp_img_%i_%i.png" %(row_nb,j)   ### HERE DRAW the substructures!!!!
                tmpfiles.append(imfile)
                worksheet.set_row(row_nb, height=size[1]) # looks like height is not in px?
                worksheet.set_column('A:B', size[0]/6.) # column width
                Draw.MolToFile(row[molCol[j]], imfile, size=size,kekulize=False)  # has to save a file on disk in order for xlsxwriter to incorporate the image
                worksheet.insert_image(row_nb, j, imfile)
        
            columns_nb = 2
            # Write data in columns and make some basic translation of numpy dtypes to excel data types
            for x in cols:
                if str(dataTypes[x]) == "object":
                    worksheet.write_string(row_nb, columns_nb, str(row[x])[:32000]) # string length is limited in xlsx
                elif ('float' in str(dataTypes[x])) or ('int' in str(dataTypes[x])):
                    if (row[x] != np.nan) or (row[x] != np.inf):
                        worksheet.write_number(row_nb, columns_nb, row[x])
                elif 'datetime' in str(dataTypes[x]):
                    worksheet.write_datetime(row_nb, columns_nb, row[x])
                columns_nb += 1
            row_nb += 1
        
        workbook.close()
        
        for f in tmpfiles:
            os.remove(f)



#############
class CalculatePvaluesContinuous:
    '''
    Calculate the p.values for each substructure in the input molecule
    '''
    def __init__(self,radii_ext):
        # needs to be the same as the one used to calcualte the dictionary of substructures
        self.radii_ext = radii_ext
        columns = ['Compound ID',
                   'Number compounds',
                   'statistic',
                   'p_value',
                   'Diff. distribution means (w - wo)',
                   'Compounds with substr.',
                   'Substructure',
                   'Substructure in Molecule'
                   ]
        self.output = pd.DataFrame(columns=columns)
        self.Bonferroni = None
    
    def calculate_p_values(self,
                           mols,
                           substructure_dictionary,
                           bioactivities,
                           mols_ids,
                           threshold_nb_substructures,
                           threshold_pvalue,
                           threshold_ratio,
                           Bonferroni=True
                           ):
        self.Bonferroni = Bonferroni
        nb_mols = float(len(set([item for sublist in substructure_dictionary.values() for item in sublist])))

        
        if type(mols) != list: mols = [mols[i] for i in np.arange(0,len(mols))]
        nb_substructures_processed=0
        already_processed = []
        for m,mol in enumerate(mols): 
            info={}
            fp = AllChem.GetMorganFingerprint(mol,max(self.radii_ext),bitInfo=info)
            substructures_sub_dict = substructure_dictionary.keys()
            
            for substructure_id, atom_radius in info.iteritems():
                nb_substructures_processed += 1
                # check if the substructure is in the reference dictionary
                if substructure_id in substructures_sub_dict and atom_radius[0][1] in self.radii_ext and substructure_id not in already_processed:
                    nb_comp_with_substructure = float(len(substructure_dictionary[substructure_id]))
                    
                    #filter threshold of compounds with the substructure
                    filter_a = nb_comp_with_substructure > threshold_nb_substructures
                    # filter threshold
                    filter_b = float(nb_comp_with_substructure / nb_mols) > threshold_ratio
                    if filter_a and filter_b:
                        mask = np.in1d(mols_ids,substructure_dictionary[substructure_id])
                        bio_substr = bioactivities[mask]
                        bio_wo_substr = bioactivities[np.logical_not(mask)]
                        # check normality
                        if sc.stats.shapiro(bio_substr) > 0.05 and sc.stats.shapiro(bio_wo_substr) > 0.05:
                            test = sc.stats.ttest_ind(bio_substr, bio_wo_substr)
                            p_value = test[1]
                            estatistic = test[0]
                        else:
                            test = sc.stats.ks_2samp(bio_substr, bio_wo_substr)
                            p_value = test[1]
                            estatistic = test[0]
                        if p_value < threshold_pvalue: 
                            env = Chem.FindAtomEnvironmentOfRadiusN(mol,atom_radius[0][1],atom_radius[0][0])
                            amap = {}
                            submol=Chem.PathToSubmol(mol,env,atomMap=amap)
                            m1 = mol
                            m1.GetSubstructMatch(submol)
                            already_processed.append(substructure_id)
                            self.output = self.output.append({
                                                             'Compound ID' : mols_ids[m],
                                                             'Number compounds':nb_mols,
                                                             'p_value': p_value,
                                                             'statistic': estatistic,
                                                             'Compounds with substr.': nb_comp_with_substructure,
                                                             'Substructure': submol,
                                                             'Substructure in Molecule': m1,
                                                             'Diff. distribution means (w - wo)': np.mean(bio_substr) - np.mean(bio_wo_substr),
                                                             },ignore_index=True)
           
        print 'number of substructures processed: ', nb_substructures_processed 
        
        if self.Bonferroni == True:
            self.output['p_value'] = self.output['p_value'] * nb_substructures_processed
            self.output = self.output[self.output.p_value < 0.05]
        
    #def HTMLOutputWriter(self,output_filename):
    #    if os.path.exists(output_filename): os.remove(output_filename)
    #    with open(output_filename, 'w') as fo:
    #        fo.write(self.output.to_html())
            
    def XlSXOutputWriter(self,frame, output_filename, molCol=['Substructure',"Substructure in Molecule"], size=(300,300)):
        """
        http://nbviewer.ipython.org/github/Team-SKI/snippets/blob/master/IPython/rdkit_hackaton/XLSX%20export.ipynb
        Saves pandas DataFrame as a xlsx file with embedded images.
        It maps numpy data types to excel cell types:
        int, float -> number
        datetime -> datetime
        object -> string (limited to 32k character - xlsx limitations)
        
        Due to xlsxwriter limitations (other python xlsx writers have the same problem) 
        temporary image files have to be written to the hard drive.
        
        Cells with compound images are a bit larger than images due to excel.
        Column width weirdness explained (from xlsxwriter docs):
        The width corresponds to the column width value that is specified in Excel. 
        It is approximately equal to the length of a string in the default font of Calibri 11. 
        Unfortunately, there is no way to specify AutoFit for a column in the Excel file format.
        This feature is only available at runtime from within Excel.
        """
        
        import xlsxwriter # don't want to make this a RDKit dependency
           
        cols = list(frame.columns)
        for i in molCol:
         cols.remove(i)

        dataTypes = dict(frame.dtypes)
        
        workbook = xlsxwriter.Workbook(output_filename) # New workbook
        worksheet = workbook.add_worksheet() # New work sheet
        worksheet.set_column('A:A', size[0]/6.) # column width
        
        # Write first row with column names
        columns_nb = 0
        for x in molCol+cols:
            worksheet.write_string(0, columns_nb, x)
            columns_nb += 1
        
        row_nb = 1
        tmpfiles = []
        for index, row in frame.iterrows():
            for j in np.arange(0,len(molCol)):
                imfile = "xlsx_tmp_img_%i_%i.png" %(row_nb,j)   ### HERE DRAW the substructures!!!!
                tmpfiles.append(imfile)
                worksheet.set_row(row_nb, height=size[1]) # looks like height is not in px?
                worksheet.set_column('A:B', size[0]/6.) # column width
                Draw.MolToFile(row[molCol[j]], imfile, size=size,kekulize=False)  # has to save a file on disk in order for xlsxwriter to incorporate the image
                worksheet.insert_image(row_nb, j, imfile)
        
            columns_nb = 2
            # Write data in columns and make some basic translation of numpy dtypes to excel data types
            for x in cols:
                if str(dataTypes[x]) == "object":
                    worksheet.write_string(row_nb, columns_nb, str(row[x])[:32000]) # string length is limited in xlsx
                elif ('float' in str(dataTypes[x])) or ('int' in str(dataTypes[x])):
                    if (row[x] != np.nan) or (row[x] != np.inf):
                        worksheet.write_number(row_nb, columns_nb, row[x])
                elif 'datetime' in str(dataTypes[x]):
                    worksheet.write_datetime(row_nb, columns_nb, row[x])
                columns_nb += 1
            row_nb += 1
        
        workbook.close()
        
        for f in tmpfiles:
            os.remove(f)
            

