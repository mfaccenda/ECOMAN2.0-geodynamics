clear

%Load the porosity distribution from your geodynamic model here
x1num = int64(?); %replace ? with the number of nodes along axis 1
x2num = int64(?); %replace ? with the number of nodes along axis 2
x3num = int64(?); %replace ? with the number of nodes along axis 3; set to 1 in 2D models

porosity = double(?); %replace ? with the 1D array containing the porosity, ordered along axis 2, axis 1, and, in 3D, axis 3.

%Write the porosity distribution
meltfilename='name of the file.h5'; %set here the porosity file name, same as in spo_input.dat
b(1)=x1num;
b(2)=x2num;
b(3)=x3num;
hdf5write(meltfilename,'Field',porosity);
hdf5write(meltfilename,'/Gridnum',b,'WriteMode','append');