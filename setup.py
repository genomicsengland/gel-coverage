from distutils.core import setup
from setuptools import find_packages

setup(
    name='gelcoverage',
    version='1.4.2',
    packages=find_packages(),
    scripts=['scripts/bigwig_analyser', 'scripts/bed_maker', 'scripts/mocked_bigwig_analyser'],
    url='',
    license='',
    author='priesgo',
    author_email='pablo.ferreiro@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'pybedtools', 'numpy', 'ujson', 'pycellbase', 'mock'],
    install_requires=[
        'pandas==0.20.3',
        'pybedtools==0.7.8',
        'ujson==1.35',
        'pyBigWig==0.3.10',
        'numpy==1.11.2',
        'pycellbase==0.3.3',
        'mock'
    ]
)
