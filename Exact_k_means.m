clear all;
clc;


x=[1,1.1,1.2,1.3,1.4,5,5.1,5.2,5.3,5.4];
N= length(x);
K=2;
mu=zeros(N,N);
eij = zeros(N,N);

for j =1:N
    for i=1:j
        eij(i,j)=0.5 .*(sum((x(i:j) - mean(x(i:j))).^2));
    end
end


% O(N^2.L) time, O(N.L) memory implementation
opplus = @opselminjoin;
idplus.val = inf;
idplus.list = {};
optimes = @opselpluscross;
idtimes.val = 0;
idtimes.list = {{}};




E = cell(N+1,1);
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

segs = E{N+1,1}.list;
    for k = 1:length(segs)
        fprintf('[');
        fprintf('(%d,%d)',segs{k}{:});
        fprintf(']');
    end
    fprintf('\n\n');
    segs = E{N+1,1}.list{1};
LE=length(segs);
E{N+1,1}.val