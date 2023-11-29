clearvars
close all

vertexs = [
    0, 0;
    5, 0;
    5, 2;
    0, 2;
    0, 0;
    ];

radius = 0.5;
theta = linspace(0,2*pi,201);
xCirc = 1 + radius*cos(theta);
yCirc = 1 + radius*sin(theta);

tempBottom = ['T = 10',176];

plot(vertexs(:,1), vertexs(:,2), '-', lineWidth = 2, color = 'black')
axis equal
axis([-0.25,5.4,-0.15,2.15])
axis off
hold on
plot(xCirc, yCirc, '-', lineWidth = 2, color = 'black')
text(-0.8,1,'$\frac{\partial T}{\partial x}\equiv 0$',...
    Interpreter='LaTex',FontSize=14)
text(5.1,1,'$\frac{\partial T}{\partial x}\equiv 0$',...
    Interpreter='LaTex',FontSize=14)
text(2.25,-0.25,'$T = 10^{\circ}$C',Interpreter='LaTeX',FontSize=12)
text(0.2,2.25,['$k_{c}\frac{\partial T}{\partial y} + \beta\left(T-T_{\infty}\right) = 0',... 
    '\quad (k_{c} = 1, \beta = 2, T_{\infty} = -5^{\circ}$C)'],...
    Interpreter='LaTeX',FontSize=12)
text(1.5,1.35,'$T = 50^{\circ}$C',Interpreter='LaTeX',FontSize=12)
hold off

print -dpng 'figureMeshPlacaForatQuad.png'