# Doublet_SW
Contains the input and python files for DARTS version 1.0.2 (TUDelft) to simulate the the doublet in a single well concept also known as recirculation well for geothermal energy production. More recent DARTS versions may require updating keywords or path of the different class and functions. More info in https://pypi.org/project/open-darts/.

Main.py: python file that structures the simulation run and post-processing. Here the model class from Model.py is called. It is used for executing DARTS.
Model.py: It is the heart of the simulation model, where the grid is created, properties are populated (Por.txt and Perm_h.txt), wells are defined, initial BC are set etc.
Por.txt: porosity array [fraction] called from Model.py.
Perm_h.txt: permeability array [mD] called from Model.py.
