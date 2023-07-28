%==========================================================================
% code for to varify the collocD.m function from the reference 185
% Greg von Winckel CGL file for find the nodes and weights
%==========================================================================
%% Generate Chebyshew Polynomial of order M=N-1 
clear;
clc;
close all;
M = 5;  % order of the polynomial
syms x;


%% Chebyshev-Gauss-Lobatto Nodes,Weights and Differential Matrix 
tau = CGL_nodes(M);  % function takes M as order of the polynomial and results LGL nodes 
nodes = tau;
D=collocD(nodes);
x_CGL = nodes
w_CGL = CGL_weights(nodes)
D_CGL = round(D,2)
D_test_CGL = CGL_Dmatrix(nodes)
% Generate the Chebyshev polynomial
cheb_poly = chebpoly(0:M, x, tau);
disp(cheb_poly);