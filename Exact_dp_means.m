clear all;
clc;
x=[1,1.1,1.2,1.3,1.4,5,5.1,5.2,5.3,5.4];
lambda=1;
N= length(x);
K=3;
mu=zeros(N,N);
eij = zeros(N,N);
eij(1,1)=lambda;
mu(1,1)=x(1);
for j=(2:N)
    for i=(1:j)
        mu(i,j) = (1/(j-i+1))*((j-i)*mu(i,j-1)+x(j));
        if i==j
             eij(i,j) = lambda;
        else
             eij(i,j)= eij(i,j-1)  + ((j-i)/(j-i+1))*(x(j)-mu(i,j-1)).^2;
        end
    end
end



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
for j = 1:N
    E{j+1} = idplus;
    for i = 1:j
        E{j+1} = opplus(E{j+1},optimes(E{i},opselinit(eij(i,j),[i,j])));
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
LE = length(segs);

E{N+1,1}.val


