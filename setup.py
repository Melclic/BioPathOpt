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
    "Bio",
    "bioservices",
    "cobra",
    "compress_json",
    "ete3",
    "networkx",
    "numpy",
    "pandas",
    "pubchempy",
    "rapidfuzz>",
    "rdkit",
    "Requests",
    "setuptools",
    "tqdm",
    "xmltodict",
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
    ], 
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
