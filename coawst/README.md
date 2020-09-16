For creating boundary, initial, and forcing files I run the master file make_coawst_files.m, which calls other matlab files, such that only make_coawst_files needs to be changed.


For making a coarse resolution grid I used [GridBuilder](https://austides.com/downloads/). I created nested grids and associated files with make_nests.m, and river forcing files with make_river_files.m.


I am running COAWST on [NeSI](https://www.nesi.org.nz). Below are some set-up notes for getting the model running on this cluster.
