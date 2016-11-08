from distutils.core import setup

setup(
    name='gelCoverage',
    version='',
    packages=['gelCoverage'],
    scripts=['coverage_summary.py'],
    url='',
    license='',
    author='mparker',
    author_email='matthew.parker@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'tqdm', 'pybedtools']
)
