from setuptools import setup
import glob

setup(
    setup_requires=["pbr>=1.8", "setuptools>=17.1"],
    pbr=True,
    scripts=glob.glob("scripts/*.py"),
    name="pyproteininference",
    version="0.9.1",
    url="",
    license="",
    author="hinklet",
    author_email="hinklet@gene.com",
    description="Python Package for running custom Protein Inference on tab-formatted tandem MS/MS search results",
)
