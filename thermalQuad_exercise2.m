clearvars
close all


tempLeft= 50.0;  %Essential B.C. for exercise 2
tempRight= 1.0;

p= [0.2, 0.8];   %Point at which temperature is interpolated in exercise 2

kc = 0.9;        %Equation coeffcients for exercise 2
a11=kc;             
a12=0.0;
a21=a12;
a22=a11;
a00=0.0;
f=0.0;

coeff=[a11,a12,a21,a22,a00,f];

meshFilesList = ['mesh4x4Quad'; 'mesh8x8Quad'];
numMeshFiles = size(meshFilesList,1);

fprintf('\nExercise 2\n')
fprintf('Interpolated temperature at point P = (%.2f, %.2f)\n',p)

for meshFileNum = 1:numMeshFiles
    meshFile = meshFilesList(meshFileNum,1:end);
    eval(meshFile);

    numNod=size(nodes,1);
    numElem=size(elem,1);

    numbering=1;
    %numbering=0;

    figure( )
    plotElementsOld(nodes,elem,numbering);
    indLeft=find(nodes(:,1) < 0.01);
    indRight=find(nodes(:,1) > 0.99);
    hold on
    plot(nodes(indLeft,1),nodes(indLeft,2),'or')
    plot(nodes(indRight,1),nodes(indRight,2),'og')
    hold off
    
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

    %Boundary Conditions (BC)
    fixedNodes=[indLeft',indRight']; %Fixed Nodes (global num.)
    freeNodes = setdiff(1:numNod,fixedNodes); %Complementary of fixed nodes

    % ------------ Natural BC
    Q=zeros(numNod,1); %initialize Q vector
    Q(freeNodes)=0;    %all of them are zero

    %------------- Essential BC
    u = zeros(numNod, 1);
    u(indLeft) = tempLeft;
    u(indRight) = tempRight;

    %------------- Reduced system
    Km=K(freeNodes,freeNodes);
    Fm=F(freeNodes) - K(freeNodes, fixedNodes)*u(fixedNodes);
    Fm=Fm+Q(freeNodes);%only for fixed nodes

    %Compute the solution
    %solve the system
    um=Km\Fm;
    u(freeNodes)=um;

    %
    %Post-process

    %Compute Q
    Q = K*u - F;

    %Plot soluition
    title = ['Solution (',meshFile(5:7),') mesh'];
    colorPalette = 'jet';
    plotContourSolution(nodes,elem,u,title,colorPalette)

    for e=1:numElem
        vertexs= nodes(elem(e,:),:);
        [alphas,isInside] = baryCoordQuad(vertexs,p);
        if (isInside > 0)
            pElem = e;
            numNodElem = elem(e,:);
            tempP = alphas*u(numNodElem);
            fprintf('For %s mesh, T(P) = %.4e\n',meshFile(5:7),tempP)
            break;
        end
    end
end 