import mogli

molecules = mogli.read('examples/dna.xyz')
mogli.BOND_RADIUS = 0.05
mogli.export(molecules[0], 'dna.html', width=1920, height=1080,
             bonds_param=1.15, camera=((40, 0, 0),
                                       (0, 0, 0),
                                       (0, 1, 0)))