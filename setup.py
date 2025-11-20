from setuptools import setup, find_packages

setup(
    name="genomextract",
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "findGenome=genomextract.findGenome:main",
            "findClosestGenome=genomextract.findClosestGenome:main",
            "assembleOrgGenome=genomextract.assembleOrgGenome:main",
            "assembleOrgGenes=genomextract.assembleOrgGenes:main",
        ]
    },
)
