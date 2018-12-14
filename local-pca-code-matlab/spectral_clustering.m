function [label1 final pick]=spectral_clustering(S1,index,K)

pick1=zeros(1,length(S1));
pick2=zeros(1,length(S1));
pick3=zeros(4,length(S1));

for i=1:length(S1)
switch class(index)
case 'double'
S2=S1{i}(index,index);
case 'cell'
S2=S1{i}(index{i},index{i});
end
X=S2;
%[U,S]=eig(S2);
%mineig=min(diag(S));
%S2=U*(S-eye(size(S1,1))*mineig)*U';
%X=U*(S-eye(size(S1,1))*mineig).^0.5;
S2_1=S2+repmat(mean(S2,1),size(S2,1),1)+repmat(mean(S2,2),1,size(S2,1))+mean(mean(S2));
rowsum=sum(S2,1);
rowsum_S2=sum(S2,1);%/(size(S2,1));
allsum_S2=sum(sum(S2));
S3=S2+1-rowsum_S2'*rowsum_S2/allsum_S2;
 
rowsum_S3=sum(S3,1);
allsum_S3=sum(sum(S3));
S4=S2+allsum_S2/size(S2,1)^2-repmat(rowsum_S2/size(S2,1),size(S2,1),1)-repmat(rowsum_S2'/size(S2,1),1,size(S2,1));
S2=diag(rowsum.^-0.5)*S2*diag(rowsum.^-0.5)+eye(size(S2,1));
rowsum=sum(S2,1);
S4=diag(rowsum.^-0.5)*S2*diag(rowsum.^-0.5);
rowsum=sum(S4,1);
S4=diag(rowsum.^-0.5)*S2*diag(rowsum.^-0.5);
%rowsum_1=sum(S2_1,1);
%S2_1=diag(rowsum_1.^-0.5)*S2_1*diag(rowsum_1.^-0.5);
%rowsum_S2_1=sum(S2_1,1)/(size(S2_1,1));
%allsum_S2_1=sum(sum(S2_1));
%S3_1=S2_1+allsum_S2_1/(size(S2_1,1))^2-repmat(rowsum_S2_1,size(S2_1,1),1)-repmat(rowsum_S2_1',1,size(S2_1,1));
 [U,S]=eig(S2);
[U,S]=svd(S2,0);

S=S-eye(size(S2));

MAXiter = 1000; % Maximum iteration for KMeans Algorithm
REPlic = 100; % Replication for KMeans Algorithm
temp=U(:,1:K)*S(1:K,1:K)^0.5;

%[label1(i,:) temp1 temp2]= kmeans(temp,K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[label1(i,:) temp1 temp2]= kmeansplusplus(temp',K,REPlic);
pick3(1,i)=sum(temp2);
%[label1(length(S1)+i,:) temp1 temp2]= kmeans(U(:,1:K),K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[label1(length(S1)+i,:) temp1 temp2]= kmeansplusplus(U(:,1:K)',K,REPlic);
pick3(2,i)=sum(temp2);
for j=1:size(temp,1)
temp(j,:)=temp(j,:)/(norm(temp(j,:))+(norm(temp(j,:))==0));
U1(j,:)=U(j,1:K)/(norm(U(j,1:K))+(norm(U(j,1:K))==0));
end

%[label1(2*length(S1)+i,:) temp1 temp2]= kmeans(temp,K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[label1(2*length(S1)+i,:) temp1 temp2]= kmeansplusplus(temp',K,REPlic);
pick3(3,i)=sum(temp2);
%[label1(3*length(S1)+i,:) temp1 temp2]= kmeans(U1,K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[label1(3*length(S1)+i,:) temp1 temp2]= kmeansplusplus(U1',K,REPlic);

pick3(4,i)=sum(temp2);
pick1(i)=S(K,K)-S(K,K);
pick2(i)=S(K,K)^0.5-S(K,K)^0.5;
%[label1(4*length(S1)+i,:)]= msgcuteig(S2,K,1);
%[U,S]=svd(S4,0);
%for j=1:size(temp,1)
%temp(j,:)=temp(j,:)/(norm(temp(j,:))+(norm(temp(j,:))==0));
%end

% [N]=size(X,1);
% U=ones(N,1)/N^0.5;
% weight=sum((X*U).^2,2).^-0.25;
% oldU=U+1;
% iter=0;
% while norm(oldU*oldU'-U*U','fro')>0.01 & iter<3
% oldU=U;
% weighted_X=diag(weight)*X*diag(weight);
% [U,S]=svd(weighted_X);%+eye(size(X,1)));

% temp=U(:,1:K)*S(1:K,1:K)*U(:,1:K)';U=U(:,1:K);
% weight=(diag(temp)./diag(weighted_X)).^-0.25;
% weight(1:5)
% iter=iter+1;

% if mod(iter,1)==0

% sum(weight.^-2)
% end
% end
% temp=U(:,1:K);
% [label1(5*length(S1)+i,:) temp1 temp2]= kmeans(temp,K,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
label1(5*length(S1)+i,:)=label1(2*length(S1)+i,:);
end
final=zeros(1,12);
[temp1 final(1)]=max(pick1);
for i=1:3
final(i+1)=i*length(S1)+final(1);
end
[temp1 final(5)]=max(pick2);
for i=1:3
final(i+5)=i*length(S1)+final(5);
end
for i=1:4
[temp1 final(8+i)]=min(pick3(i,:));
final(8+i)=final(8+i)+(i-1)*length(S1);
end
pick=[max(pick1),max(pick2),min([pick3]')];

