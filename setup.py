from distutils.core import setup

setup(
    name='gelcoverage',
    version='1.1.1',
    packages=['gelcoverage'],
    scripts=['scripts/bigwig_analyser'],
    url='',
    license='',
    author='priesgo',
    author_email='pablo.ferreiro@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'pybedtools', 'numpy', 'ujson']
)
