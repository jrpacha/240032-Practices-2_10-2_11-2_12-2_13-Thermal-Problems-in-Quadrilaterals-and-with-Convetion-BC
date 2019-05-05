function [W,pt2D]=gaussValues2DQuad(n)
% 2D from the 1D values
[wx,ptx]=gaussValues1D(n);
[wy,pty]=gaussValues1D(n);
%points and weights 
pt2D=[];
W=[];
for i=1:n
    for j=1:n
        pt2D=[pt2D; [ptx(i),pty(j)]];
        W=[W,wx(i)*wy(j)];
    end
end
