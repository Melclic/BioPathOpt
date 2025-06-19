import os
from setuptools import setup, find_packages

major_version = 0
minor_version = 0
patch_version = 1

if os.environ.get("PATCH_NUMBER"):
    patch_version = os.environ.get("PATCH_NUMBER")

version = str(major_version) + "." + str(minor_version) + "." + str(patch_version)

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, "README.md")).read()

third_party_require = [
    "Bio>=1.8.0",
    "bioservices>=1.12.1",
    "cobra>=0.29.1",
    "compress_json>=1.1.1",
    "ete3>=3.1.3",
    "networkx>=3.4.2",
    "numpy>=1.26.4",
    "pandas>=2.3.0",
    "pubchempy>=1.0.4",
    "rapidfuzz>=3.13.0",
    "rdkit>=2025.3.3",
    "Requests>=2.32.4",
    "scikit_learn>=1.6.1",
    "setuptools>=69.5.1",
    "torch>=2.5.1",
    "tqdm>=4.66.4",
    "xmltodict>=0.14.2",
]

tests_require = [
    "nosexcover",
    "coverage",
    "nose-timer",
    "nose-xunitmp",
    "pylint",
    "nose",
]

require = third_party_require + tests_require


setup(
    name="biopathopt",
    version=version,
    description="BioPathOpt is a package to optimize and explore metabolic engineering",
    long_description=README,
    classifiers=[
        "Programming Language :: Python :: 3.10",
    ],  # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords="",
    author="Melchior du Lac",
    author_email="",
    license="",
    packages=find_packages(exclude=["ez_setup"]),
    package_dir={},
    package_data={},
    include_package_data=True,
    zip_safe=False,
    test_suite="nose.collector",
    install_requires=require,
    tests_require=tests_require,
)
