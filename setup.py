from setuptools import setup
import glob

setup(
    setup_requires=["pbr>=1.8", "setuptools>=17.1"],
    pbr=True,
    scripts=glob.glob("scripts/*.py"),
    name="pyproteininference",
    version="0.9.8",
    url="",
    license="",
    author="hinklet",
    author_email="hinklet@gene.com",
    description="Python Package for running protein inference algorithms on tab-formatted tandem MS/MS search results.",
    keywords=['protein inference', 'proteomics', 'mass spectrometry'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    long_description=long_description,
    long_description_content_type='text/markdown',
)
