%% Engine Cooling fin
clearvars
close all

kc= 0.5;         %Thermal conductivity
ff= 0.0;
tempRight= 90.0; %Temperature at the right boundary
beta= 0.05;      %Convection at top, bottom and left boundaries
Tinf=-17.0; 
p=[0.9,2.1];     %Point to interpolate the temperature      

%Load geometry: node coordinates and elements
eval('meshAleta2DQuad')
numNod= size(nodes,1);
numElem= size(elem,1);

numbering=0;
plotElements(nodes, elem, numbering);

indTop= find(nodes(:,2)-2*nodes(:,1)/5 > 1.99); %top boundary's nodes' indices
indLeft=find(nodes(:,1) < 0.01); %left boundary's nodes' indices
indRight=find(nodes(:,1) > 4.99) %right boundary's nodes' indices
%Id for the circle's nodes

hold on
plot(nodes(indTop,1),nodes(indTop,2),'o','color','black',...
    'markerFaceColor','green','markerSize',5)
plot(nodes(indLeft,1),nodes(indLeft,2),'o','color','black',...
    'markerFaceColor','green','markerSize',5)
plot(nodes(indRight,1),nodes(indRight,2),'o','color','black',...
    'markerFaceColor','red','markerSize',5)
hold off

%Define Coefficients vector of the model equation
%In this case we use the Poisson coefficients defined in the problem above
a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=ff;
coeff=[a11,a12,a21,a22,a00,f];

K=zeros(numNod);
F=zeros(numNod,1);
Q=zeros(numNod,1);

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
end %end for elements
%we save a copy of the initial F array
%for the postprocess step
Kini= K;
Fini= F;

%Booundary Conditions
fixedNodes= indRight'; %fixed Nodes (global numbering)
freeNodes= setdiff(1:numNod,fixedNodes); %free Nodes (global numbering)

%------------- Convetion BC
indCV=[indTop',indLeft']; %must be a row vector
indCV=unique(indCV);      %just to avoid duplicate indices
[K,Q]=applyConvQuad(indCV,beta,Tinf,K,Q,nodes,elem);

% ------------ Essential BC
u=zeros(numNod,1); %initialize uu vector
u(indRight)=tempRight;
Fm=F(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);%here u can be 
                                                      %different from zero 
                                                      %only for fixed nodes
%Reduced system
Km=K(freeNodes,freeNodes);
Fm=Fm+Q(freeNodes);

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
% format long e
% ts=table(int16(tableSol(:,1)),tableSol(:,2),tableSol(:,3),tableSol(:,4),...
%     tableSol(:,5),'variableNames',{'NumNod','X','Y','T','Q'});
% writetable(ts,fileName);

%Compute the temperature at point p=[0.9,2.1].

for e=1:numElem
    vertexs= nodes(elem(e,:),:);
    [alphas,isInside] = baryCoordQuad(vertexs,p);
    if (isInside > 0)
        pElem = e;
        numNodElem= elem(e,:);
        break;
    end
end

tempP = alphas*u(numNodElem);
fprintf('\nExercise 2:\n')
fprintf('Point P = (%.1f,%.1f) belongs to element number: %d\n',p,pElem)
fprintf('Number of nodes of elem %d: %d, %d, %d, %d\n',pElem,numNodElem)
fprintf('Interpolated temperature at point P: %.4e\n',tempP)
