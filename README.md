# mogli

mogli is a python module for easily rendering molecules.

The fastest way to get started is using mogli.py in the command line:
```
% python -m mogli <file_name> 
```

This will open a new window in which the molecule will be rendered. You can then rotate the molecule using your mouse.

In python, you can achieve something similar with three lines of code:

```
>>> import mogli
>>> molecules = mogli.read('examples/dna.xyz')
>>> mogli.show(molecules[0])
```

![DNA rendered in a GR window](https://raw.githubusercontent.com/FlorianRhiem/mogli/doc-images/dna-cli.png)

## Atomic bonds
You might notice that the atomic bonds between the atoms don't look right. You can hide the bonds by adding the `show_bonds=False` to your `show()` call, but of course there's a a better way to fix the missing bonds here.

For the example above, adding `bonds_param=1.15` to `show()` will do the trick. mogli currently offers two ways of calculating atomic bonds. The first method compares the distance between every two atoms with the sum of their valence radii. If they are further apart, no bond is formed. To allow a bit of adjustment, this radii sum is multiplied with a factor that can be set using the `bonds_param` parameter to `show()` and `draw()`. By default, `1.0` is used.

The second method uses a constant maximum distance instead. If you use this method by passing `bonds_method='constant_delta'` to `show()` or `draw()`, you can set the constant distance with the `bonds_param` parameter.

## Exporting to files
Instead of viewing molecules interactively, you can export the molecule as well, for example as Portable Network Graphics into a `.png` file, or as HTML5 page with interactive WebGL code as `.html` file. To do so, simply call the `export()` function, like so:

```
>>> import mogli
>>> molecules = mogli.read('examples/dna.xyz')
>>> mogli.export(molecules[0], 'dna.html', width=1920, height=1080,
                 bonds_param=1.15, camera=((40, 0, 0),
                                           (0, 0, 0),
                                           (0, 1, 0)))
```

This example code also shows setting the camera in code by passing the `camera` parameter. It expects three sequences of three floating point numbers: the camera position, the center of focus and a direction that should point upward. *Exporting to HTML5 with WebGL is a pretty near feature, give it a try!*

## Drawing with GR
In case you use the GR framework, you can use mogli to draw molecules into your GR graphics. To do so, just call `draw()`. You can use the parameters `xmin`, `xmax`, `ymin`, `ymax`, `width` and `height` just like you would when using `gr.drawimage()`.

```
>>> import gr
>>> import mogli
>>> molecules = mogli.read('examples/dna.xyz')
>>> mogli.draw(molecules[0])
>>> gr.updatews()
```

![DNA rendered in a GR window](https://raw.githubusercontent.com/FlorianRhiem/mogli/doc-images/dna-gr.png)

## Dependencies
mogli depends on GR3, which is included in the [GR framework](http://gr-framework.org/) ([PyPI]( https://pypi.python.org/pypi/gr), [GitHub](https://github.com/jheinen/gr)), on [glfw3](http://www.glfw.org/) ([GitHub](https://github.com/glfw/glfw)) and python bindings for glfw3 ([PyPI](https://pypi.python.org/pypi/glfw), [GitHub](https://github.com/FlorianRhiem/pyGLFW)).

Additionally, [Pybel](http://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html) should be installed; mogli will work without it, but only for xyz files.
