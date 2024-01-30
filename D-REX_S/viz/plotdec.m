function plotdec()

global printmod input_dir output_dir figname strain

v = [0      0.4470 0.7410;
     0.8500 0.3250 0.0980;
     0.9290 0.6940 0.1250;
     0.4940 0.1840 0.5560;
     0.4660 0.6740 0.1880;
     0.3010 0.7450 0.9330;
     0.6350 0.0780 0.1840];

filename = [input_dir,'anisdec.h5'];

perc_anis = h5read(filename,'/perc_anis');
perc_hexa = h5read(filename,'/perc_hexa');
perc_tetra= h5read(filename,'/perc_tetra');
perc_ortho= h5read(filename,'/perc_ortho');
perc_mono = h5read(filename,'/perc_mono');
perc_tri  = h5read(filename,'/perc_tri');

f = figure('Name','Tensor decomposition','NumberTitle','off'); clf;
subplot(211)
plot(strain,perc_anis,'-','Color',v(1,:),'LineWidth',1.5);
title('Norm of total anisotropy / Norm full tensor  (%)')
xlabel('\gamma strain ')
ylabel('(%)')

subplot(212)
plot(strain,perc_hexa,'-','Color',v(2,:),'LineWidth',1.5);
hold on
plot(strain,perc_tetra,'-','Color',v(3,:),'LineWidth',1.5);
plot(strain,perc_ortho,'-','Color',v(4,:),'LineWidth',1.5);
plot(strain,perc_mono,'-','Color',v(5,:),'LineWidth',1.5);
plot(strain,perc_tri,'-','Color',v(6,:),'LineWidth',1.5);
hold off
title('Norm of * anisotropy / Norm total anisotropy (%)')
xlabel('\gamma strain ')
ylabel('%')
legend('hexagonal', 'tetragonal', 'orthorhombic','monoclinic','triclinic','Location','NorthWest')

if printmod
    print('-dpng', '-r150',[output_dir,figname]);
end