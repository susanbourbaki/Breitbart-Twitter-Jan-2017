% calcGraphLaplacian -- function that computes the un-normalized graph
% Laplacian.
%
% INPUTS:
%
% A -- adjacency matrix, with A(i,j) weight of arc from j to i.
%
% OUTPUTS:
% L -- unnormalized Laplacian, L = W-A.
% w -- vector of out-degrees.  W = diag(w).
%
function [L,w] = calcGraphLaplacian(A)

w = sum(A,1);  % column sums
L = diag(w) - A;


