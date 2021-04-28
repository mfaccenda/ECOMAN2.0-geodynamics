%Make vtp*.h5 files for cartesian or polar/spherical chunk domains
%In 3D, ignore the 3rd dimensions
clear

%Number of nodes
nx1= 1  ;%number of nodes along axis 1
nx2= 1  ;%number of nodes along axis 2
nx3= 1  ;%number of nodes along axis 3, set to 1 in 2D
nodenum=nx1*nx2*nx3;


%File infos
finit = 1; %Initial file number
fstp  = 1; %Step file number
fend  = 1; %Final file number
for f = finit:fstp:fend
    
    %A) Load here the Eulerian fields of your geodynamic model
    V1model=...; velocity along the axis 1
    V2model=...; velocity along the axis 2
    V3model=...; velocity along the axis 3
    Tkmodel=...; Temperature
    Pmodel=... ; Pressure
    Fdmodel=...; Fraction of deformation accommodated by disl. creep (adimensional)
        
    time=...; %total elapsed time
    dt=...  ; %elapsed time from previous output file
        
    %B) Convert your arrays into 1D D-REX_M format
    %for example, if your array are in 2 or 3 dimensions
    
    V1=zeros(nodenum,1); V2=V1; V3=V1; Tk=V1; P=V1; Fd=V1;
    for i1 = 1:nx1
        for i2 = 1:nx2
            for i3 = 1:nx3
                
                %Global index for D-REX_M 1D array
                i4 = i2 + (i1-1)*nx2 + (i3-1)*nx1*nx2;
                
                V1(i4) = V1model(i1,i2,i3); % in m/s or adimensional
                V2(i4) = V2model(i1,i2,i3); % in m/s or adimensional
                V3(i4) = V3model(i1,i2,i3); % in m/s or adimensional
                Tk(i4) = Tkmodel(i1,i2,i3); % in Kelvin
                P(i4)  =  Pmodel(i1,i2,i3); % in Pa
                Fd(i4) = Fdmodel(i1,i2,i3); % adimensional, 0?Fd?1
            end
        end
    end
    
    t=[dt time]; % array with time infos (in sec or adimensional)
    
    %C)Save arrays in HDF5 format files
    fname=['vtp', num2str(f,'%.4d'),'.h5'];
    hdf5write(fname,'/Nodes/V1',V1);
    hdf5write(fname,'/Nodes/V2',V2, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/V3',V3, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/Tk',Tk, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/P',P, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/Fd',Fd, 'WriteMode', 'append');
    h5writeatt(fname,'/','Time',t);
    
end

