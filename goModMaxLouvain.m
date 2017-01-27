% goModMaxLouvain -- script which maximizes modularity using Louvain
% algorithm.
%
% INPUTS:
%
% A -- weighted adjacency matrix.  Convention: A(i,j) = weight of arc from j to i
% names -- names associated with different vertices.
%
% OUTPUTS:
% 
% S -- identified communities.
% Q -- modularity.

function [S,Q] = goModMaxLouvain(A,names)

w_out = sum(A,1);  % column sums of A
w_in = sum(A,2);   % row sums of A
w = sum(w_out);    % sum of all arc weights

B = A - w_out'*w_in'/w;

% Louvain to maximize modularity
[S,Q] = genlouvain(B);

ncomm = max(S);
for i=1:ncomm
    inds = find(S == i);
    disp(['Community ',int2str(i),': ']);
    for j=1:length(inds)
        disp([names{inds(j)}]);
    end
    disp('*******************');
    %pause;
end
