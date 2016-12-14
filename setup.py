from distutils.core import setup

setup(
    name='gelcoverage',
    version='1.1.0',
    packages=['gelcoverage'],
    scripts=['scripts/bigwig_analyser.py'],
    url='',
    license='',
    author='priesgo',
    author_email='pablo.ferreiro@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'pybedtools', 'numpy', 'ujson']
)
