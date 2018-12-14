function [label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(X,M,d,size_of_neighborhood,neighborhood_type)
%input: X data set. M: correct label. 
%d: dimension. 
%size_of_neighborhood: the size of neighborhood in the covering 
%neighborhood_type: different ways of generating covering, can be set as 1.1, which is the one used in the simulations of the paper "Spectral Clustering Based on Local PCA".
%When neighborhood_type = 5, it means fixed K_NN neighborhood for every point. When neighborhood_type = 6, it means neighborhood with fixed radius size for every point.
% label: output label. allerror: misclassification rate
K=max(M);
[ nets Dimss]=covering_map(X,size_of_neighborhood,neighborhood_type);
for i=1:length(nets)
models{i}=find_model(X(nets{i},:),1,Dimss);
end
for i=1:length(nets)
[U,S]=svd(models{i}.cov);
models{i}.cov3=U(:,1:d)*U(:,1:d)';
end
for i=1:length(nets)
net_truce(i)=M(nets{i}(1));
end
for i=1:length(nets)
for j=1:length(nets)
A(i,j)=norm(models{i}.mean-models{j}.mean);
B(i,j)=norm(models{i}.cov3-models{j}.cov3);
%B(i,j)=norm(models{i}.cov-models{j}.cov);
end
end
% B=B./A;
% for i=1:size(B,1)
% B(i,i)=0;
% end
% C=zeros(size(A));
% for i=1:length(nets)
% [bb,bb1]=sort(A(i,:),'ascend');
% C(i,bb1(1:10))=1;
% end 
% C=((C+C')>0.5);
C=ones(size(A));
minA=min(min(A+10^3*eye(size(A))));
maxA=max(max(A+0*eye(size(A))));
minB=min(min(B+10^3*eye(size(B))));
maxB=max(max(B+0*eye(size(B))));
stepA=(log(maxA)-log(minA))/10;
sigmaA=exp([-10:1:10]*stepA+log(minA));
stepB=(log(maxB)-log(minB))/10;
sigmaB=exp([-2:0.4:10]*stepB+log(minB));
%%%two circles
%if minA==0
connected_A=max(min(A+100*eye(size(A,1))));
sigmaA=2.^([0:0.5:0])*connected_A;
connected_graph=(A<connected_A);
connected_graph=connected_graph-diag(diag(connected_graph));
if sum(sum(connected_graph))>0
temp=reshape(B.*connected_graph,1,size(B,1)^2);
average_eta=quantile(temp(find(temp>0)),0.5);

else
average_eta=1;
end
disp(average_eta)
disp(sigmaA)
%sigmaA=exp([-7:1:3])*maxA;
%else
%stepA=(log(maxA)-log(minA))/10;
%sigmaA=exp([-10:1:10]*stepA+log(minA));
%end
%if minB==0
sigmaB=exp([-4:0.5:0]);%*maxB;
%else
%stepB=(log(maxB)-log(minB))/10;
%sigmaB=exp([-7:0.4:10]*stepB+log(minB));
%end
allerror1=ones(length(sigmaA),length(sigmaB),4);
besterror=inf;
for i=1:length(sigmaA)
sigmaB=2.^([0:0.5:0])*average_eta;
allSigmaB(i,:)=sigmaB;
for j=1:length(sigmaB)
S1{1}=exp(-A.^2/sigmaA(i)^2).*exp(-B.^2/sigmaB(j)^2);
S1{1}=S1{1}.*C; 
disp(S1)
cd ./spectral' clustering'/
[label1]=SpectralClustering(S1{1},K,1);
cd ..
for k=1:4
allerror1(i,j,k)=testerror(label1(k,:),net_truce);
end

if allerror1(i,j,4)<besterror
min_i=i;min_j=j;
bestlabel=label1;
besterror=allerror1(i,j,4);
end
end

end
label1=bestlabel;
% temp=allerror1(:,:,4);
% [n1,n2]=size(temp);
% [aa,aa1]=min(reshape(temp,1,n1*n2));
% min_i=mod(aa1-1,n1)+1;
% min_j=floor((aa1-1)/n1)+1;
% i=min_i;j=min_j;
% S1{1}=exp(-A.^2/sigmaA(i)^2).*exp(-B.^2/sigmaB(j)^2);
% S1{1}=S1{1}.*C; 
% cd ../spectral' clustering'/
% [label1]=SpectralClustering(S1{1},K,1);
% cd ../newproject
[mapping]=nearest_model(X,models,0,1,0);%find the distances between the points and the models
label=zeros(size(label1,1),size(X,1));
for i=1:size(label1,1)
label(i,:)=find_label(X,label1(i,:),mapping);%assign each point to the nearest model
end

allerror=zeros(1,size(label1,1));
for i=1:length(allerror)
allerror(i)=testerror(label(i,:),M)/length(M);
end

coarse_label=zeros(size(label1,1),size(X,1));
for j=1:size(label1,1)
for i=1:length(nets)
coarse_label(j,nets{i})=label1(j,i);
end
end

disp(average_eta)
