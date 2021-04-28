clear;

clf

path='/Users/manuelefaccenda/Desktop/DREX/08052020/Fallingslab/';
fname=[path,'syntomo0020.h5'];

X1=h5read(fname,'/X1'); nx1=size(X1,1);
X2=h5read(fname,'/X2'); nx2=size(X2,1);
X3=h5read(fname,'/X3'); nx3=size(X3,1);
nodenum = nx1*nx2*nx3;
Rho=h5read(fname,'/Rho');
X1save=h5read(fname,'/Sav');

Sav=zeros(nodenum,6,6);
m=0;
for i1=1:nx1
    for i2=1:nx2
        for i3=1:nx3
            
            gi = i1 + (i2-1)*nx1 + (i3-1)*nx1*nx2;
            m = m + 1;
            Sav(gi,1,1) = X1save(1,m);
            Sav(gi,2,2) = X1save(2,m);
            Sav(gi,3,3) = X1save(3,m);
            Sav(gi,2,3) = X1save(4,m);
            Sav(gi,3,2) = X1save(4,m);
            Sav(gi,1,3) = X1save(5,m);
            Sav(gi,3,1) = X1save(5,m);
            Sav(gi,1,2) = X1save(6,m);
            Sav(gi,2,1) = X1save(6,m);
            Sav(gi,4,4) = X1save(7,m);
            Sav(gi,5,5) = X1save(8,m);
            Sav(gi,6,6) = X1save(9,m);
            Sav(gi,1,4) = X1save(10,m);
            Sav(gi,4,1) = X1save(10,m);
            Sav(gi,2,5) = X1save(11,m);
            Sav(gi,5,2) = X1save(11,m);
            Sav(gi,3,6) = X1save(12,m);
            Sav(gi,6,3) = X1save(12,m);
            Sav(gi,3,4) = X1save(13,m);
            Sav(gi,4,3) = X1save(13,m);
            Sav(gi,1,5) = X1save(14,m);
            Sav(gi,5,1) = X1save(14,m);
            Sav(gi,2,6) = X1save(15,m);
            Sav(gi,6,2) = X1save(15,m);
            Sav(gi,2,4) = X1save(16,m);
            Sav(gi,4,2) = X1save(16,m);
            Sav(gi,3,5) = X1save(17,m);
            Sav(gi,5,3) = X1save(17,m);
            Sav(gi,1,6) = X1save(18,m);
            Sav(gi,6,1) = X1save(18,m);
            Sav(gi,5,6) = X1save(19,m);
            Sav(gi,6,5) = X1save(19,m);
            Sav(gi,4,6) = X1save(20,m);
            Sav(gi,6,4) = X1save(20,m);
            Sav(gi,4,5) = X1save(21,m);
            Sav(gi,5,4) = X1save(21,m);
            
        end
    end
end
