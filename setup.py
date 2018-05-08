from setuptools import setup, find_packages

setup(
    name="modpin",  
   	version="1.1.0", 
	author="GRIB",
	author_email="", 
    description="SET OF PYTHON SCRIPTS TO MODEL AND ANALYZE PROTEIN-PROTEIN INTERACTIONS",
	keywords="protein interaction modelling",
	url="",
    classifiers=[  
        "Development Status :: 3 - Alpha",
    ],

    packages=find_packages(),
    include_package_data = True,
	install_requires=[
   		'numpy',
   		'biopython',
		'scipy',
		'matplotlib'
	],
 	
	entry_points = {
        'console_scripts': ['configuration=modpin.config_setup:main',
'modppi=modpin.scripts.modppi:main',
'analysis=modpin.scripts.analysis:main',
'modelist=modpin.scripts.modelist:main',
'FilterHomologPairs=modpin.scripts.FilterHomologPairs:main',
'Renumerate=modpin.scripts.Renumerate:main',
'3did = modpin.src.functions.3did:main',

],
	}

	

)





