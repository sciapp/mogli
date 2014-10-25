import gr
import mogli

molecules = mogli.read('examples/dna.xyz')
mogli.ATOM_RADIUS = 3
gr.clearws()
gr.setviewport(0, 1, 0, 1)
mogli.draw(molecules[0], bonds_param=1.15, camera=((60, 0, 0),
                                                    (0, 0, 0),
                                                    (0, 1, 0)))
gr.settextalign(gr.TEXT_HALIGN_CENTER, gr.TEXT_VALIGN_TOP)
gr.text(0.5, 1.0, 'DNA example')
gr.updatews()
