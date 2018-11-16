function [label1] = SpectralClustering(S2,K,opt)
rowsum=sum(S2,1);
S2=diag(rowsum.^-0.5)*S2*diag(rowsum.^-0.5)+eye(size(S2,1));
[U,S]=svd(S2,0);
S=S-eye(size(S2));
MAXiter = 1000; % Maximum iteration for KMeans Algorithm
REPlic = 1000; % Replication for KMeans Algorithm
temp=U(:,1:K)*S(1:K,1:K)^0.5;
if opt==1
[label1(1,:) temp1 temp2]= kmeans(temp,K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[label1(2,:) temp1 temp2]= kmeans(U(:,1:K),K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
else
[label1(1,:) temp1 temp2]= kmeansplusplus(temp',K,REPlic);
[label1(2,:) temp1 temp2]= kmeansplusplus(U(:,1:K)',K,REPlic);
end

for j=1:size(temp,1)
temp(j,:)=temp(j,:)/(norm(temp(j,:))+(norm(temp(j,:))==0));
U1(j,:)=U(j,1:K)/(norm(U(j,1:K))+(norm(U(j,1:K))==0));
end
if opt==1
[label1(3,:) temp1 temp2]= kmeans(temp,K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[label1(4,:) temp1 temp2]= kmeans(U1,K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
else
    [label1(3,:) temp1 temp2]= kmeansplusplus(temp',K,REPlic);
    [label1(4,:) temp1 temp2]= kmeansplusplus(U1',K,REPlic);
end