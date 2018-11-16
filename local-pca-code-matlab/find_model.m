function model=find_model(X,model_opt,Dimss)
model.mean=mean(X,1);
%model.mean=X(1,:);
model.num=size(X,1);
%model.cov=cov(X-repmat(model.mean,size(X,1),1));%+eye(size(X,2));
Y=X-repmat(model.mean,size(X,1),1);
switch model_opt
case 1
model.cov=Y'*Y/(size(Y,1)-1);
[U,S]=svd(model.cov,0);
diagS=diag(S)';
%if length(diagS)<size(X,2)
%diagS((length(diagS)+1):size(X,2))=0;
%end
model.cov=model.cov;%/trace(model.cov);%+50*eye(size(X,2));
model.cov2=model.cov;%+0.000001*eye(size(X,2));
model.dim=size(X,2);
case 1.1
model.cov=Y'*Y/(size(Y,1)-1);
[U,S]=svd(model.cov,0);
diagS=diag(S)';
%if length(diagS)<size(X,2)
%diagS((length(diagS)+1):size(X,2))=0;
%end
model.cov=U(:,1:2)*U(:,1:2)';%+50*eye(size(X,2));
model.cov2=model.cov;%+0.000001*eye(size(X,2));
model.dim=size(X,2);
case 2 
model.cov=Y'*Y/(size(Y,1)-1);
[U,S]=svd(model.cov,0);
diagS=diag(S)';
if length(diagS)<size(X,2)
diagS((length(diagS)+1):size(X,2))=0;
end
%model.cov2=U(:,find(diagS>((mean(diagS))/2)));
model.cov2=U(:,1:Dimss);
model.dim=max(size(model.cov2,2),1);
model.cov2=model.cov2*model.cov2';%+0.0001*eye(size(X,2));
%model.cov2=U(:,model.dim)*S(model.dim,model.dim)*U(:,model.dim)';
 model.cov=model.cov2;
case 3
Z=zeros(size(X,1),size(X,2)*(size(X,2)+3)/2);
Z(:,1:size(X,2))=X;
pointer=size(X,2)+1;
for i=1:size(X,2)
for j=1:i
Z(:,pointer)=X(:,i).*X(:,j);
pointer=pointer+1;
end
end
model.mean=mean(Z,1);
model.num=size(X,1);
Y=Z-repmat(model.mean,model.num,1);

model.cov=Y'*Y/(size(Y,1)-1);
[U,S]=svd(model.cov,0);
diagS=diag(S)';
model.cov2=U(:,find(diagS>((mean(diagS))/2)));
model.dim=size(model.cov2,2);
model.cov2=model.cov2*model.cov2';
model.cov=model.cov2;
case 5
model.cov=eye(size(X,2));

model.cov2=model.cov;%+0.000001*eye(size(X,2));
model.dim=size(X,2);
case 6
model.cov=Y'*Y/(size(Y,1)-1);
[U,S]=svd(model.cov,0);
diagS=diag(S)';
if length(diagS)<size(X,2)
diagS((length(diagS)+1):size(X,2))=0;
end
model.cov2=U(:,find(diagS>((mean(diagS))/2)));
model.dim=max(size(model.cov2,2),1);
model.noise=(diagS(min(model.dim+1,length(diagS))))^0.5;
model.cov2=model.cov;%+0.000001*eye(size(X,2));
model.dim=size(X,2);
%model.proj=U(:,find(diagS>((mean(diagS))/2)));
model.proj=U(:,1:2);
model.dim=size(model.proj,2);
model.proj=model.proj*model.proj';
model.cov2=model.proj;
end  
