% calcNetworkRisk -- function which computes the 'network risk' associated
% with nodes in a strongly connected network.
%
% INPUTS:
%
% L -- unnormalized graph Laplacian.
%
% OUTPUTS:
% 
% u -- normalized basis for nullspace of L, such that u>0 and sum(u)=1.

function u = calcNetworkRisk(L)

u = null(L);
if u(1)<0
    u = -u;
end

u = u./sum(u);