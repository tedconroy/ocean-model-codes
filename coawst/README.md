For creating boundary, initial, and forcing files I run the master file make_coawst_files.m, which calls other matlab files, such that only make_coawst_files needs to be changed.


For making a coarse resolution grid I used [GridBuilder](https://austides.com/downloads/). I created nested grids and associated files with make_nests.m, and river forcing files with make_river_files.m.


I am running COAWST on [NeSI](https://www.nesi.org.nz). Below are some set-up notes for getting the model running on this cluster. First obtain the source code (can use svn as well).

```
git clone https://github.com/jcwarner-usgs/COAWST
```
in my bash.rc file I put the module info needed for COAWST
```
module load netCDF-Fortran/4.5.2-gimpi-2020a
```
To make SCRIP, I set FORT = gfortran in the makefile, and had to change the lib and include directories to somewhere that I had permission to write. I then set the env vars to 
```
export   MCT_INCDIR=/home/tc196/COAWST/Lib/MCT/include 
export   MCT_LIBDIR=/home/tc196/COAWST/Lib/MCT/lib 
```
to find my Netcdf env, i used commands nc-config -all, which nc-config, nc-config --flibs, to find
```
export  NETCDF_INCDIR=/opt/nesi/CS400_centos7_bdw/netCDF-Fortran/4.5.2-gimpi-2020a/include 
export  NETCDF_LIBDIR=/opt/nesi/CS400_centos7_bdw/netCDF-Fortran/4.5.2-gimpi-2020a/lib 
export  NETCDF=/opt/nesi/CS400_centos7_bdw/netCDF-Fortran/4.5.2-gimpi-2020a 
export  NETCDF_CONFIG=/opt/nesi/CS400_centos7_bdw/netCDF-Fortran/4.5.2-gimpi-2020a/bin /nc-config 
 ```
 now update the coawst.bash file accordingly, and compile COAWST
 ```
 ./coawst.bash
 ```

