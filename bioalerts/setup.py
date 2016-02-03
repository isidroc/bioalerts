"""
Set up the bioalerts package
Authors: Isidro Cortes-Ciriano <isidrolauscher@gmail.com>
Date: 27.10.14
"""

from distutils.core import setup 

#packagedata={'bioalerts': ['FP_Calculator.py','BioAlerts.py']}
#packagedata={'propy': ['aaindex1','aaindex2','aaindex3','html/*','instruction/*','data/*','aaindex/*']}


setup(name = 'bioalerts', 
	version = '0.1.0', 
	description ="  ",
	author = "Isidro Cortes-Ciriano",
	author_email = "isidrolauscher@gmail.com",
	license = "GPL",
	#packages = ['bioalerts'],
	#package_data=packagedata,
#	data_files = datafiles,
	#package_dir={''},
	scripts = ['Alerts.py','FPCalculator.py','LoadMolecules.py'], 
#	py_modules = [],
#	setup_requires=["numpy"],
	install_requires = ['numpy','scipy','operator','pandas']
)

