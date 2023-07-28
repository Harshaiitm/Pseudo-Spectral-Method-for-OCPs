function [z,D]=legDc(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% legDc.m
%
% Computes the Legendre differentiation matrix with collocation at the 
% Legendre-Gauss-Lobatto nodes.
%
% Reference: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 05/26/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncation + 1
N1=N+1;
% CGL nodes
zc=cos(pi*(0:N)/N)';
% Uniform nodes
zu=linspace(-1,1,N1)';
% Make a close first guess to reduce iterations
if N<3
    z=zc;
else
    z=zc+sin(pi*zu)./(4*N);
end
% The Legendre Vandermonde Matrix
P=zeros(N1,N1);
% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.
zold=2;
while max(abs(z-zold))>eps
    zold=z;
        
    P(:,1)=1;    P(:,2)=z;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*z.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    z=zold-( z.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
end
Z=repmat(z,1,N1);
Zdiff=Z-Z'+eye(N1);
L=repmat(P(:,N1),1,N1);
L(1:(N1+1):N1*N1)=1;
D=(L./(Zdiff.*L'));
D(1:(N1+1):N1*N1)=0;
D(1)=(N1*N)/4;
D(N1*N1)=-(N1*N)/4;   
