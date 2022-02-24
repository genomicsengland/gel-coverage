# from distutils.core import setup
from setuptools import find_packages, setup
import os
import io
import re

# read the contents of your README file

VERSION = '1.4.6'

REGEX_COMMENT = re.compile(r"[\s^]#(.*)")


def parse_requirements(filename):
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    with open(filename, "rt") as filehandle:
        return tuple(
            filter(None, (REGEX_COMMENT.sub("", line).strip() for line in filehandle))
        )

setup(
    name='gel-coverage',
    version=VERSION,
    packages=find_packages(),
    scripts=['scripts/bigwig_analyser', 'scripts/bed_maker'],
    url='https://github.com/genomicsengland/gel-coverage',
    download_url="https://github.com/genomicsengland/gel-coverage/archive/v{}.tar.gz".format(VERSION),
    license='Apache',
    author='Luca Venturini',
    author_email='luca.venturini@genomicsengland.co.uk',
    description='Whole genome coverage analysis tool',
    install_requires=parse_requirements("requirements.txt"),
    extras_require={"test": ["httpretty>=0.9.5", "unittest"]},
    classifiers=[
        'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7'
      ]
)
