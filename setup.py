from distutils.core import setup
from setuptools import find_packages

setup(
    name='gelcoverage',
    version='1.2.3',
    packages=find_packages(),
    scripts=['scripts/bigwig_analyser'],
    url='',
    license='',
    author='priesgo',
    author_email='pablo.ferreiro@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'pybedtools', 'numpy', 'ujson', 'pycellbase'],
    install_requires=[
        'pandas',
        'pybedtools==0.7.8',
        'ujson==1.33',
        'pyBigWig==0.3.4',
        'numpy==1.10.4',
        'pycellbase==0.3.2'
    ]
)
