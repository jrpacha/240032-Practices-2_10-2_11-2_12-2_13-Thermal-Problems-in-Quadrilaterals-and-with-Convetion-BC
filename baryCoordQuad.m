function [alphas, isInside]=baryCoordQuad(vertices, point)
% (c)Numerical Factory 2021
isInside=1; 
v1=vertices(1,:);
v2=vertices(2,:);
v3=vertices(3,:);
v4=vertices(4,:);
% define the previous notation from the interpolation formula:
% point=(1-mu)*(1-lambda)*v1+lambda*(1-mu)*v2+mu*lambda*v3+(1-lambda)*mu*v4;
% arranging terms:
% 0 = a + b*lambda + c*mu + d*lambda*mu; 
a=v1-point;
b=v2-v1;
c=v4-v1;
d=v1-v2-v4+v3;

% we use (x,y) instead of (lambda, mu) for easy code
% solution is supposed to have 0 <= x <=1 and 0 <= y <=1 
x=0.5; y=0.5; %initial point
dx=0.1*ones(2,1);
iter=0;
tol=1.e-12;
while(norm(dx) > tol && iter < 20 ) 
    %apply Newton-Rapson method to solve f(x,y)=0
    f = a + b*x + c*y + d*x*y; 
    % Newton: x_{n+1} = x_n - (Df^-1)*f
    % or equivalently denoting dx = x_{n+1}-x_n
    % Newton: Df*dx=-f
    Df(:,1)=(b + d*y);  %df/dx
    Df(:,2)=(c + d*x);  %df/dy
    bb=-f'; %independent term
    dx=Df\bb;
    x=x+dx(1);
    y=y+dx(2);
    iter=iter+1;
    if (norm([x,y]) > 10), iter=20; end %non convergent: just to save time
end
if (iter < 20)
    lambda=x;
    mu=y;
    %
    alpha1=(1-mu)*(1-lambda);
    alpha2=lambda*(1-mu);
    alpha3=mu*lambda;
    alpha4=(1-lambda)*mu;
    alphas=[alpha1,alpha2,alpha3,alpha4];
    if ( min(alphas) < 0)
        isInside=0;
    end
else
    alphas=[-1,-1,-1,-1];
    isInside=0;
end