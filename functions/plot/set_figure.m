function f = set_figure(m,n)
    f = figure;
    myPaperSize = [m n];
    lowerLeft = [0 5];
    set(f,'PaperUnits','centimeters');
    set(f,'PaperPositionMode','manual');
    set(f,'PaperSize',myPaperSize);
    set(f,'Units','centimeters');
    set(f,'Position',[lowerLeft myPaperSize]);
    set(f,'PaperPosition',[lowerLeft myPaperSize]);
%     
%     set(f,'defaultAxesFontSize',8)
%     set(f, 'defaultTextInterpreter','latex');
%     set(f, 'defaultAxesTickLabelInterpreter','latex');
%     set(f, 'defaultLegendInterpreter','latex');
%     
%     set(f, 'defaultLegendItemTokenSize',[5,5]);
%     set(groot, 'defaultLegendItemTokenSize',[10,1]);

end