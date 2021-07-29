import re
from setuptools import setup, find_packages
def get_version():
	try:
		f = open("ssDripPipelines/_version.py")
	except EnvironmentError:
		return None
	for line in f:
		line=line.rstrip()
		m = re.match("__version__ = '([^']+)'", line)
		if m:
			ver = m.group(1)
			return ver
	return None

setup(
	name='ssDripPipeline',
	version=get_version(),
	author='LiKuan',
	author_email="396777306@qq.com",
	packages=find_packages(),
	scripts=['bin/ssDRIPSeqAnalysis.py',],
	#include_package_data=True,
	url='https://github.com/PEHGP/ssDripPipeline',
	license='GPL-3.0',
	description='Useful tools for the analysis of ssDRIP-seq data ',
	long_description='',
	classifiers=[
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics'],
	install_requires=[
		"pandas >= 1.0.1",
		"biopython >=1.78",
	],
	zip_safe=True,
)
