from distutils.core import setup

setup(
    name='gelcoverage',
    version='1.0.0',
    packages=['gelcoverage'],
    scripts=['scripts/coverage_summary'],
    url='',
    license='',
    author='mparker',
    author_email='matthew.parker@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'tqdm', 'pybedtools']
)
