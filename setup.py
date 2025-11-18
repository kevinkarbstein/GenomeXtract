from setuptools import setup

setup(
    name="ncbi-genome-tools",
    version="0.1.0",
    packages=["felix"],
    entry_points={
        "console_scripts": [
            "findGenome = felix.findGenome:main",
            "findClosestGenome = felix.findClosestGenome:main",
            "assembleGenome = felix.assembleGenome:main",
            "assembleGenes = felix.assembleGenes:main",
        ]
    },
)

