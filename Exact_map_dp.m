clear all;
clc;
x=[1,1.1,1.2,1.3,1.4,5,5.1,5.2,5.3,5.4];
N = length(x);

alpha0=0.05;
sigma0=(1)*std(x);
sigma_hat=(0.4)*std(x);
m0=mean(x);

for j = 1:N
    for i = 1:j
        sigma_k=(1/sigma0^2+length(x(i:j))/sigma_hat^2)^(-1);
        mu_k=sigma_k*(m0/sigma0^2+sum(x(i:j))/sigma_hat^2);
        eij(i,j) = 0;
        for q=i:j
            eij(i,j)= eij(i,j) + Gaussian_nll(x(q),mu_k,sigma_k,sigma_hat);
        end
        eij(i,j) = eij(i,j) - log(alpha0) - gammaln(length(x(i:j)));
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

E{N+1,1}.val-gammaln(alpha0) + gammaln(alpha0+N)
 
 
 
function nl = Gaussian_nll(X,mu_k,sigma_k,sigma_hat)
   nl = (1/(2*(sigma_k+sigma_hat^2)))*sum((abs(X-mu_k)).^2)+(1/2)*log(sigma_k+sigma_hat^2); 
end

