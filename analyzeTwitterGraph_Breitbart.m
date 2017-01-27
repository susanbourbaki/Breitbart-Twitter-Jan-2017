% analyzeTwitterGraph_Breitbart -- analysis for Jan 11 blog post on
% susanbourbaki.com

addpath ./GenLouvain2;
load 'twitter_Breitbart.mat';


%% Community detection via modularity maximization and Louvain

% News sources only
disp('-----------------------------------------');
disp('Community detection -- news sources only');
disp('-----------------------------------------');
[Sn,Qn] = goModMaxLouvain(A(news_inds,news_inds),seeds.names(news_inds));

% Full graph
disp('-----------------------------------------');
disp('Community detection -- all nodes');
disp('-----------------------------------------');
[S,Q] = goModMaxLouvain(A,seeds.names);

%% Centrality within the extremist community
% 'degree centrality' -- number of followers
[ds,di] = sort(seeds.nfollowers(extremist_inds),'descend');
ee_degcent = seeds.names(extremist_inds(di));

Aee = A(extremist_inds,extremist_inds);

% PageRank without teleportation 
% displays Rank (value of the stationary distribution)
Le = calcGraphLaplacian(Aee);  
u = calcNetworkRisk(Le);
[us,iu] = sort(u,'descend');
ee_pagerank = seeds.names(extremist_inds(iu));
disp('--------------------------------');
disp('Centrality for white nationalist nodes.');
disp('PageRank without teleportation:');
disp('--------------------------------');
for i=1:length(iu)
    disp([seeds.names{extremist_inds(iu(i))},' -- ',int2str(i),' (',num2str(us(i)),')']);
end

%% News preferences for white nationalists
Ane = A(news_inds,extremist_inds);
[mval,minds] = max(Ane);
prefnews = seeds.names(news_inds(minds));  % 6 prefer Breitbart; 3 prefer Drudge

%% Number of unique followers of white nationalist focal nodes
%% List of all unique followers of white nationalist accounts

n = length(extremist_inds);

fol_wn = [];
nfol = 0;

for i=1:n
    curid = seeds.inds(extremist_inds(i));
    foli = find(edges(:,1) == curid);
    
    nfol = nfol + length(foli);
    fol_wn = [fol_wn; edges(foli,2)];
end

fol_wn = unique(fol_wn);
n_fol_wn = length(fol_wn);

disp(['Sum of all followers: ',int2str(nfol)]);
disp(['Number of unique followers of extremist nodes: ',int2str(n_fol_wn)]);


%% Co-followership between white nationalists and Trump

trump_ind = 23;
trump_id = seeds.inds(trump_ind);

trump_foli = find(edges(:,1) == trump_id);
trump_fol = edges(trump_foli,2);

trump_cofol = intersect(trump_fol,fol_wn);

disp(['Number of Trump followers who follow any extremist node: ',int2str(length(trump_cofol))]);
prop_WNforTrump = length(trump_cofol)/n_fol_wn;
str = ['Proportion of white nationalist followers who also follow Trump:'];
disp([str,num2str(prop_WNforTrump)]);

% Weighted average of arc weights pointing from white nationalists to
% Trump, weighted by PageRank:
wn_avg_prop = A(trump_ind,extremist_inds)*u;


%% Calculate news overlaps -- e.g. followers of white nationalists who also follow both Breitbart and Drudge, etc
%
% n -- number of white nationalist focal nodes
%
% Matrix: n by 3+3+1
% columns:
%   1 -- Breitbart
%   2 -- Drudge
%   3 -- Fox
%   4 -- Breitbart + Drudge
%   5 -- Breitbart + Fox
%   6 -- Drudge + Fox
%   7 -- Breitbart + Drudge + Fox

news_subinds = [1,9,12];  % news_inds(news_subinds) = Breitbart, Drudge, Fox
news_ids = seeds.inds(news_inds(news_subinds));

breit_foli = find(edges(:,1) == news_ids(1));
drudge_foli = find(edges(:,1) == news_ids(2));
fox_foli = find(edges(:,1) == news_ids(3));

breit_fol = edges(breit_foli,2);
drudge_fol = edges(drudge_foli,2);
fox_fol = edges(fox_foli,2);

n = length(extremist_inds);
% columns:
%   1 -- Breitbart only
%   2 -- Drudge only
%   3 -- Fox only
%   4 -- Breitbart + Drudge, but not Fox
%   5 -- Breitbart + Fox, but not Drudge
%   6 -- Drudge + Fox, but not Breitbart
%   7 -- Breitbart + Drudge + Fox

overlaps = zeros(n,7);

for i=1:n
%   % find indices of followers
    disp(['Current node: ',seeds.names{extremist_inds(i)}]);

    curid = seeds.inds(extremist_inds(i));
    foli = find(edges(:,1) == curid);
    fol = edges(foli,2);
    
    disp(['Total followers: ',int2str(length(fol))]);
    
    % all mutual followers with Breitbart
    ball = intersect(fol,breit_fol);
    disp(['Followers who also follow Breitbart: ',int2str(length(ball))]);

    % all mutual followers with Drudge
    dall = intersect(fol,drudge_fol);
    disp(['Followers who also follow Drudge: ',int2str(length(dall))]);

    % all mutual followers with Fox
    fall = intersect(fol,fox_fol);
    disp(['Followers who also follow Fox: ',int2str(length(fall))]);
    
    % all mutual followers between Breitbart and Drudge
    bdall = intersect(ball,dall);
        
    % all mutual followers between Breitbart and Fox
    bfall = intersect(ball,fall);
    
    % all mutual followers between Drudge and Fox
    dfall = intersect(dall,fall);
    
    % Followers of all three
    bdf = intersect(bdall,fall);    
    
    nall = length(bdf);
    
    nbd = length(bdall) - nall;
    nbf = length(bfall) - nall;
    ndf = length(dfall) - nall;
    
    nb = length(ball) - nbd - nbf - nall;
    nd = length(dall) - nbd - ndf - nall;
    nf = length(fall) - nbf - ndf - nall;
    
    overlaps(i,1) = nb;
    overlaps(i,2) = nd;
    overlaps(i,3) = nf;
    overlaps(i,4) = nbd;
    overlaps(i,5) = nbf;
    overlaps(i,6) = ndf;
    overlaps(i,7) = nall;
    
 end
