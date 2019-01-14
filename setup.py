import os
from setuptools import setup
import sys
import versioneer


if sys.version_info < (3, 6):
    sys.stdout.write("ngsindex requires python3.6.\n")
    sys.exit(1)


with open(
    os.path.join(os.path.abspath(os.path.dirname(__file__)), "README.md"),
    encoding="utf-8"
) as f:
    long_description = f.read()


setup(
    name='ngsindex',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Utilities for working with NGS index formats.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/jdidion/ngsindex',
    author='John Didion',
    author_email='github@didion.net',
    license='MIT',
    packages=['ngsindex'],
    install_requires=[
        'xphyle>=4.0.0rc0'
    ],
    tests_require=[
        'pytest',
        'pytest-cov',
        'pytest-datadir',
        'dataclasses'
    ],
    entry_points={
        "console_script": [
            "summarize-index=ngsindex.cli:summarize"
        ]
    },
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
