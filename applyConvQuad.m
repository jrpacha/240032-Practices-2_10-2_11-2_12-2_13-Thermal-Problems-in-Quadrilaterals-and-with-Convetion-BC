function [K,Q]=applyConvQuad(indCV,beta,Tinf,K,Q,nodes,elem)
%
% (c) Numerical Factory 2018
%
%------------------------------------------------------------------------
% Apply a convection BC on a boundary of a quadrilateral meshed domain. 
%
% Input:
% indCV -> List of nodes on the convection boundary. They can be
%          non-connected but it needs: 
%          i) Each connected component contains at least two nodes 
%          ii) Only one corner of a quadrilateral is accepted (two edges of the same quadrilateral). 
%          If an unique quadrilateral is part of the CV bondary, we have two call twice 
%          this function: one time with a corner and then the rest, and sume contributions
%
%   *      *           *                *          *              *
%   |      |           |                |          |              |
%   |      |     =     |          +     |   o be   |    +         |
%   *------*           *------*       ..*          *..     *------*
%
%
%  Hint: See the Triang version for the rest of variables and more explanations
%------------------------------------------------------------------------------

numElem=size(elem,1); 
numCo=length(indCV);
if (numCo==1) 
    error('applyConvQuad: It can''t manage a single node'); 
end
Kaux=(beta/6)*[2 1; 1 2]; 
Faux=0.5*beta*Tinf*[1; 1]; 
for k=1:numElem
  aux=[0,0,0,0];
  for j=1:4
    r=find(indCV==elem(k,j)); %find if j node of element k is in the convection nodes list 
      if(~isempty(r)), aux(j)=1; end
  end
  %See the Triang version for the explanation of the binary number    
  number=aux*[1;2;4;8];
  % Here we codify both edges and corners on the boundaries
  switch (number) 
    case  3, ij=[1,2,0,0];  % edge 1:   aux=[1,1,0,0];
    case  6, ij=[2,3,0,0];  % edge 2:   aux=[0,1,1,0];
    case  7, ij=[1,2,2,3];  % corner 2: aux=[1,1,1,0];
    case  9, ij=[4,1,0,0];  % edge 4
    case 11, ij=[4,1,1,2];  % corner 1 
    case 12, ij=[3,4,0,0];  % edge 3
    case 13, ij=[3,4,4,1];  % corner 4
    case 14, ij=[2,3,3,4];  % corner 3
    case 15, error('applyConvQuad: It can''t manage two corners'); % two corners 
    otherwise, ij=[0,0,0,0];
  end
  for j=1:2:3
   if ~ij(j)==0
    n1=elem(k,ij(j)); 
    n2=elem(k,ij(j+1));
    h=norm(nodes(n1,:)-nodes(n2,:)); 
    fico=[n1, n2];   
    K(fico,fico)=K(fico,fico)+h*Kaux; 
    Q(fico)=Q(fico)+h*Faux;
   end
  end
end

