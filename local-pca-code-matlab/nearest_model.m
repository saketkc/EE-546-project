function [mapping]=nearest_model(X,models,regularization,model_type,nets,eps)
mapping=zeros(size(X,1),length(models));
switch model_type
case 1
for i=1:length(models)
mapping(:,i)=sum(models{i}.mean.^2)+sum(X.^2,2)-2*X*models{i}.mean';
end
case 1.1
Z=zeros(size(X,1),size(X,2)*(size(X,2)+3)/2);
Z(:,1:size(X,2))=X;
pointer=size(X,2)+1;
for i=1:size(X,2)
for j=1:i
Z(:,pointer)=X(:,i).*X(:,j);
pointer=pointer+1;
end
end
for i=1:length(models)
if regularization(i)==0 & models{i}.dim<size(models{i}.cov2)
if norm((models{i}.mean-models{j}.mean)*(eye(size(models{i}.cov2,1))-models{i}.cov2))==0
S1(i,j)=(models{i}.mean-models{j}.mean)*(models{i}.mean-models{j}.mean)';
else
S1(i,j)=((models{i}.mean-models{j}.mean)*(models{i}.cov2+regularization(i)*eye(size(models{i}.cov2,1)))^(-1)*(models{i}.mean-models{j}.mean)');
end
else
for j=1:size(X,1)

mapping(j,i)=((Z(j,:)-models{i}.mean)*(models{i}.cov+regularization(i)*eye(size(models{i}.cov,1)))^(-1)*(Z(j,:)-models{i}.mean)');
end
end
end
otherwise
%regu=zeros(1,length(models));

%for j=1:length(models)
%clear dist
%dist=zeros(1,length(nets{j}));
%for jj=1:length(nets{j})
%dist(jj)=((X(nets{j}(jj),:)-models{j}.mean)*(-models{j}.cov+eye(size(models{j}.cov,1)))*(X(nets{j}(jj),:)-models{j}.mean)');
%end
%regu(j)=max(real(median(dist.^0.5)^2/(eps^2)),0);
%end
%regu
for j=1:length(models)
if abs(regularization(j))<10^-5 & models{j}.dim<size(models{j}.cov,1)
mapping(:,j)=max(sum(((X-repmat(models{j}.mean,size(X,1),1))*(models{j}.cov+regularization(j)*eye(size(models{j}.cov,1)))^(-0.5)).^2,2).^0.5,sum(((X-repmat(models{j}.mean,size(X,1),1))*(eye(size(models{j}.cov,1)))^(-0.5)).^2,2).^0.5);

else
mapping(:,j)=sum(((X-repmat(models{j}.mean,size(X,1),1))*(models{j}.cov+regularization(j)*eye(size(models{j}.cov,1)))^(-0.5)).^2,2).^0.5;
end
end
end