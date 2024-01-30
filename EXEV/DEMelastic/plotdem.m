function plotdem(x,y,z,a1a2num,a2a3num,Title,logscale)

for i=1:a1a2num
    z(i,i+1:end)=z(i,i);
end

pcolor(x,y,z)
colorbar
shading interp
axis square
if logscale
    %     pcolor(x,y,log10(z))
    % colorbar
    % shading interp
    % axis square
    set(gca, 'XScale', 'log', 'YScale', 'log');
    %set(gca,'YMinorTick','on');
    %set(gca,'Layer','top');
end
grid on
title(Title,'fontsize',20,'fontweight','bold','fontangle','italic')

hold on
if logscale
    [C,h1]=contour(x,y,z,'k','LineWidth',1.0);
    xlabel('a2/a3','fontsize',14,'fontweight','bold','fontangle','italic')
    ylabel('a1/a2','fontsize',14,'fontweight','bold','fontangle','italic')
else
    [C,h1]=contour(log10(x),log10(y),z,'k','LineWidth',1.0);
    xlabel('log(a2/a3)','fontsize',14,'fontweight','bold','fontangle','italic')
    ylabel('log(a1/a2)','fontsize',14,'fontweight','bold','fontangle','italic')

end

clabel(C,h1);
hold off

h=gca;
set(h,'fontsize',12,'fontweight','bold','fontangle','italic')
c = colorbar;%('southoutside');
%c.Label.String = Title;
set(c,'fontsize',12,'fontweight','bold','fontangle','italic')


hold on
xr=(0:1:10000);
plot(xr,xr,'r','LineWidth',1.5)
hold off
