from distutils.core import setup

setup(
    name='gelCoverage',
    version='1.0.0',
    packages=['gelCoverage'],
    scripts=['scripts/coverage_summary'],
    url='',
    license='',
    author='mparker',
    author_email='matthew.parker@genomicsengland.co.uk',
    description='', requires=['pandas', 'pyBigWig', 'tqdm', 'pybedtools']
)
