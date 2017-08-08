function SetFigureDefaults(w,h)
% Figure size, font and axes placement
% Usage:
%       SetFigureDefaults(width,height)
%       width and height in cm
%
sfX=0.75;sfY=0.75;
set(0,'DefaultAxesPosition',[0.15,0.15,sfX,sfY])
set(0,'DefaultAxesFontName','Times')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultTextFontSize',6)

set(gcf,'PaperUnits','centimeters','Units','Centimeters')
p1=get(gcf,'Position');
p2=get(gcf,'PaperPosition');
p1(2)=p1(2)+p1(4)-h/sfY;
p1([3,4])=[w/sfX,h/sfY];
p2([3,4])=[w/sfX,h/sfY];
set(gcf,'Position',p1,'PaperPosition',p2)