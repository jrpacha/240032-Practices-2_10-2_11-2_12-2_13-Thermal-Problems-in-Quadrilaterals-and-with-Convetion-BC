clearvars
close all

fileName='solutionThermalQuad.xlsx';

tempRight = 0.0;
tempLeft = 0.0;

eval('mesh2x2Quad');

numNod=size(nodes,1);
numElem=size(elem,1);
numbering=1;
%numbering=0;
%figure()
plotElementsOld(nodes,elem,numbering);

indLeft=find(nodes(:,1) < 0.01);
indRight=find(nodes(:,1) > 0.99);
hold on
plot(nodes(indLeft,1),nodes(indLeft,2),'or')
plot(nodes(indRight,1),nodes(indRight,2),'og')
hold off

a11=1;
a12=1;
a21=a12;
%a21=2;
a22=a11;
a00=1;
f=1;

coeff=[a11,a12,a21,a22,a00,f];
K=zeros(numNod);
F=zeros(numNod,1);
for e=1:numElem
    [Ke,Fe]=bilinearQuadElement(coeff,nodes,elem,e);
    %
    % Assemble the elements
    %
    rows=[elem(e,1); elem(e,2); elem(e,3); elem(e,4)];
    colums= rows;
    K(rows,colums)=K(rows,colums)+Ke; %assembly
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end 
Kini = K;
Fini = F;

%Boundary Conditions (BC)
%Impose Boundary Conditions for this example. In this case only essential 
%and natural BC are considered. We will do a thermal example for convection 
%BC (mixed).
%indRigth=[3;6;9];
%indLeft=[1;4;7];

fixedNodes=[indLeft',indRight']; %Fixed Nodes (global num.)
freeNodes = setdiff(1:numNod,fixedNodes); %Complementary of fixed nodes

% ------------ Natural BC
Q=zeros(numNod,1); %initialize Q vector
Q(freeNodes)=0; %all of them are zero


% ------------ Essential BC
u=zeros(numNod,1); %initialize u vector
u(indRight)=tempRight; %all of them are zero
u(indLeft)=tempLeft;
Fm=F(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);%here u can be 
                                                      %different from zero 
                                                      % modify the linear system, this is only valid if BC = 0.
%Reduced system
Km=K(freeNodes,freeNodes);
Fm=Fm+Q(freeNodes);%only for free nodes
%Compute the solution
%solve the System


format short e; %just to a better view of the numbers
um=Km\Fm;
u(freeNodes)=um;
%u

%PostProcess: Compute secondary variables and plot results
Q=Kini*u-Fini;
titol='Equation solution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale)

%Fancy output
tableSol=[(1:numNod)',nodes,u,Q];
fprintf('%8s%9s%15s%15s%14s\n','Num.Nod','X','Y','T','Q')
fprintf('%5d%18.7e%15.7e%15.7e%15.7e\n',tableSol')

%write an Excel with the solutions
format long e
ts=table(int16(tableSol(:,1)),tableSol(:,2),tableSol(:,3),tableSol(:,4),...
     tableSol(:,5),'variableNames',{'NumNod','X','Y','T','Q'});
writetable(ts,fileName);