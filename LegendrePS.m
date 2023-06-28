%==========================================================================
% code for to varify the collocD.m function from the reference 185
% Greg von Winckel LGL,LGR,LG file for find the nodes and weights
%==========================================================================
%% Generate Legendre Polynomial of order M=N-1 
clear;
clc;
close all;
N = 6;  % Number of nodes
syms x;
Legendre_Polynomial =  legendreP(N,x)

%% Legendre-Gauss-Lobatto Nodes,Weights and Differential Matrix 
[nodes,weights] = LGL_nodes(N);  % function takes N as number of nodes and results LGL nodes 
[x,Dmat]=legDc(N);
D=collocD(nodes);
x_LGL = round(nodes,2)
w_LGL = weights
D_LGL = round(D,2)
D_test_LGL = round(Dmat,2)

%% Legendre-Gauss-Radau Nodes,Weights and Differential Matrix 
[nodes,weights] = LGR_nodes(N);
x_LGR = round(nodes,2)
w_LGR = weights
D=collocD(nodes);
D_LGR = round(D,2)

%% Legendre-Gauss Quadrature Weights and Nodes
a = -1;            % interval (-1,1) 
b = 1;
[nodes,weights] = LG_nodes(N,a,b);
x_LG = round(nodes,2)
w_LG = weights;
D=collocD(nodes);
D_LG = round(D,2)

% plot the Polynomial
% hold on;
% for M = 0:5      % order of the polynomial M=N-1
%     Legendre_Polynomial = legendreP(M, x_LG);  % Evaluate Legendre polynomial of degree N for x_LG
%     plot(x_LG,Legendre_Polynomial,'DisplayName', ['M = ' num2str(M)]);  % Plot the Legendre polynomial
% end
% 
% hold off;
% legend('show');  % Show legend with the polynomial degree
% xlabel('x_LG');  % Label the x-axis
% ylabel('Legendre Polynomial');  % Label the y-axis
% title('Legendre Polynomials');  % Add a title to the plot

