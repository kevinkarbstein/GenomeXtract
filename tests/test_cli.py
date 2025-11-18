import subprocess

def test_findGenome_help():
    subprocess.run(["findGenome", "--help"], check=True)

def test_findClosestGenome_help():
    subprocess.run(["findClosestGenome", "--help"], check=True)

def test_assembleOrgGenome_help():
    subprocess.run(["assembleOrgGenome", "--help"], check=True)

def test_assembleOrgGenes_help():
    subprocess.run(["assembleOrgGenes", "--help"], check=True)
