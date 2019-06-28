# from distutils.core import setup
from setuptools import find_packages, setup
import os


# read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

VERSION = '1.4.2'

setup(
    name='gelcoverage',
    version=VERSION,
    packages=find_packages(),
    scripts=['scripts/bigwig_analyser', 'scripts/bed_maker'],
    url='https://github.com/genomicsengland/gel-coverage',
    download_url="https://github.com/genomicsengland/gel-coverage/archive/v{}.tar.gz".format(VERSION),
    license='Apache',
    author='Pablo Riesgo Ferreiro',
    author_email='pablo.riesgo-ferreiro@genomicsengland.co.uk',
    description='Whole genome coverage analysis tool',
    long_description=long_description,
    long_description_content_type='text/markdown',
    requires=['pandas', 'pyBigWig', 'pybedtools', 'numpy', 'ujson', 'pycellbase'],
    install_requires=[
        'pandas==0.20.3',
        'pybedtools==0.7.8',
        'ujson==1.35',
        'pyBigWig==0.3.4',
        'numpy==1.11.2',
        'pycellbase==0.3.3'
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7'
      ]
)
