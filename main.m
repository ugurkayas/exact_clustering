clc;
clear;

x=[5.3,1,5.1,5.4,1.2,1.4,1.3,5,5.2,1.1];           % sample data

k=2;                                               % the number of clusters

lambda=1;                                          % penalizer for dp-means      

alpha0=0.05;                                       % hyper-parameters for MAP-DP
sigma0=(1)*std(x);
sigma_hat=(0.5)*std(x);
m0=mean(x);

[y1,c1]= Exact_k_means(x,k);

[y2,c2]= Exact_dp_means(x,lambda);

[y3,c3]= Exact_map_dp(x,alpha0,sigma0,sigma_hat,m0);

% Display results

fprintf('K-means:\n');
fprintf('Error Value: %f\n', y1);
fprintf('Cluster Setup: ');
disp(c1);

fprintf('DP-means:\n');
fprintf('Error Value: %f\n', y2);
fprintf('Cluster Setup: ');
disp(c2);

fprintf('MAP-DP:\n');
fprintf('Error Value: %f\n', y3);
fprintf('Cluster Setup: ');
disp(c3);
