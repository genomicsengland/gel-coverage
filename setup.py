from distutils.core import setup
from setuptools import find_packages

setup(
    name='gelcoverage',
    version='1.1.2',
    packages=find_packages(),
    scripts=['scripts/bigwig_analyser'],
    url='',
    license='',
    author='priesgo',
    author_email='pablo.ferreiro@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'pybedtools', 'numpy', 'ujson']
)
