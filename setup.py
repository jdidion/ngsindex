from setuptools import setup
import sys
import versioneer


version_info = sys.version_info
if version_info < (3, 6):
    sys.stdout.write("ngsindex requires python3.6.\n")
    sys.exit(1)


setup(
    name='ngsindex',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Utilities for working with NGS index formats.',
    url='https://github.com/jdidion/ngsindex',
    author='John Didion',
    author_email='github@didion.net',
    license='MIT',
    packages=['ngsindex'],
    install_requires=[
        'xphyle',
        'pysam'
    ],
    tests_require=[
        'pytest',
        'pytest-cov',
        'pytest-datadir',
        'dataclasses'
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'License :: OSI Approved :: MIT License',
        'License :: Public Domain',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
    ],
)
