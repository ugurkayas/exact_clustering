clear all;
clc;


x=[1,1.1,1.2,1.3,1.4,5,5.1,5.2,5.3,5.4];                                              % sample data
K=2                                                                                   % the number of clusters
N= length(x);                                                                         % length of the data
mu=zeros(N,N);                                                                        % global mean matris for each possible mean from i to j
eij = zeros(N,N);                                                                     % objective value matris for each possible mean from i to j

for j =1:N                                                                            % objective value matris calculation
    for i=1:j
        eij(i,j)=0.5 .*(sum((x(i:j) - mean(x(i:j))).^2));
    end
end


% O(N^2.L) time, O(N.L) memory implementation                                         % min-plus semiring with tuples
opplus = @opselminjoin;
idplus.val = inf;
idplus.list = {};
optimes = @opselpluscross;
idtimes.val = 0;
idtimes.list = {{}};




E = cell(N+1,1);                                                                         % Lifted Bellman's recursion
for i = 1:N+1
    E{i,1} = idplus;
end
E{1,1} = idtimes;
for m=1:K                                                                              
  for j = N:-1:1
   for i = 1:j
        E{j+1} = opplus(E{j+1},optimes(E{i},opselinit(eij(i,j),[i,j])));
    end
  end
end

segs = E{N+1,1}.list;                                                                  % extracting the final clustering setting
for k = 1:length(segs)
    fprintf('[');
    fprintf('(%d,%d)',segs{k}{:});
    fprintf(']');
end
fprintf('\n\n');                                                                       % the final clustering setting
segs = E{N+1,1}.list{1};
LE = length(segs);

E{N+1,1}.val                                                                           % objective value

