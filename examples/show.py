import mogli

molecules = mogli.read('examples/dna.xyz')
for molecule in molecules:
    mogli.show(molecule, bonds_param=1.15)