function [x,w,P]=LGR_nodes(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LGR_nodes.m
%
% Computes the Legendre-Gauss-Radau nodes, weights and the LGR Vandermonde 
% matrix. The LGR nodes are the zeros of P_N(x)+P_{N+1}(x). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncation + 1
N1=N+1;
% Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes
x=-cos(2*pi*(0:N)/(2*N+1))';
% The Legendre Vandermonde Matrix
P=zeros(N1,N1+1);
% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.
xold=2;
% Free abscissae
free=2:N1;
while max(abs(x-xold))>eps
    xold=x;
     
    P(1,:)=(-1).^(0:N1);       
    
    P(free,1)=1;    P(free,2)=x(free);
       
    for k=2:N1
        P(free,k+1)=( (2*k-1)*x(free).*P(free,k)-(k-1)*P(free,k-1) )/k;
    end
     
    x(free)=xold(free)-((1-xold(free))/N1).*(P(free,N1)+P(free,N1+1))...
                ./(P(free,N1)-P(free,N1+1));
        end
nodes =double(x(free));
% The Legendre-Gauss-Radau Vandermonde
P=P(1:N1,1:N1);
% Compute the weights
w=zeros(N1,1);
w(1)=2/N1^2;
w(free)=(1-x(free))./(N1*P(free,N1)).^2;