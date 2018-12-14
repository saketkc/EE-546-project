function [S1 neighbor]=find_distances(models,D,metric,regularization,sigma,eps,nets,X,size_of_neighborhood)
n=length(models);
S=zeros(n,n);
S1=zeros(n,n);
size_of_neighborhood=min(size_of_neighborhood,n);
%size_of_neighborhood=floor(1.5*log(n));
sigma=min(sigma,n);

for i=1:n
    if abs(regularization(i))<10^-5 & models{i}.dim<size(models{i}.cov,1)
        
        for j=1:n
            S1(i,j)=max(((models{i}.mean-models{j}.mean)*(models{i}.cov2+regularization(i)*eye(size(models{i}.cov2,1)))^(-1)*(models{i}.mean-models{j}.mean)'),(models{i}.mean-models{j}.mean)*(models{i}.mean-models{j}.mean)');
        end
    else
        for j=1:n
            %S(i,j)=norm(models{j}.mean-models{i}.mean);
            S1(i,j)=((models{i}.mean-models{j}.mean)*(models{i}.cov2+regularization(i)*eye(size(models{i}.cov2,1)))^(-1)*(models{i}.mean-models{j}.mean)');
        end
    end
end
%S1=(S1+S1')/2;
size_of_neighborhood=10; 
for i=1:n
    [temp1 temp]=sort(S1(i,:));
    neighbor{i}=temp(1:size_of_neighborhood);
end


S3=zeros(n,n);
S4=zeros(n,n);

%for i=1:n
%for j=1:n
%    %S(i,j)=norm(models{j}.mean-models{i}.mean);
%    S3(i,j)=((models{i}.mean-models{j}.mean)*(models{i}.cov2)*(models{i}.mean-models{j}.mean)');
%	S4(i,j)=((models{i}.mean-models{j}.mean)*(eye(size(models{i}.cov2,1))-models{i}.cov2)*(models{i}.mean-models{j}.mean)');
%end
%end


%for i=1:n
%clear dist
%for jj=1:length(nets{i})
%dist(jj)=((X(nets{i}(jj),:)-models{i}.mean)*(-models{i}.cov2+eye(size(models{i}.cov2,1)))*(X(nets{i}(jj),:)-models{i}.mean)');
%end
%S5(i)=median(dist);
%length(find(S4(i,:)>50*median(dist)))
%S3(i,find(S4(i,:)>50*median(dist)))=inf;
%[temp1 temp]=sort(S3(i,:));
%neighbor{i}=temp(1:size_of_neighborhood);
%end

for i=1:n
    for j=neighbor{i}
        switch metric
            case 1
                S(i,j)=norm((models{j}.cov2*logm(models{j}.cov2^-0.5*models{i}.cov2*models{j}.cov2^-0.5)*models{j}.cov2^0.5),'fro');
            case 2
                S(i,j)=norm(models{j}.cov2-models{i}.cov2,'fro');
            case 3% hellinger
                S(i,j)=(1-2^(D/2)*det(models{i}.cov2)^(1/4)*det(models{j}.cov2)^(1/4)/det(models{i}.cov2+models{j}.cov2)^0.5)^0.5;
            case 4
                
                S(i,j)=(2^(D/2)*det(models{i}.cov2+regularization*eye(D))^(1/4)*det(models{j}.cov2+regularization*eye(D))^(1/4)/det(models{i}.cov2+models{j}.cov2+2*regularization*eye(D))^0.5)*exp(-0.25*(models{i}.mean-models{j}.mean)*(models{i}.cov2+models{j}.cov2+2*regularization*eye(D))^-1*(models{i}.mean-models{j}.mean)');
            case 5
                S(i,j)=norm(models{i}.cov2*(models{i}.num-1)-models{j}.cov2*(models{j}.num-1),'fro')^2/((norm(models{i}.cov2*(models{i}.num-1),'fro')^2)+(norm(models{j}.cov2*(models{j}.num-1),'fro')^2));
            case 6%Hellinger
                S(i,j)=(2^(D/2)*det(models{i}.cov2)^(1/4)*det(models{j}.cov2)^(1/4)/det(models{i}.cov2+models{j}.cov2)^0.5)*exp(-0.25*(models{i}.mean-models{j}.mean)*(models{i}.cov2+models{j}.cov2)^-1*(models{i}.mean-models{j}.mean)');
            case 6.5
			S(i,j)=exp(-norm(models{i}.mean-models{j}.mean)^2/eps^2)*exp(-norm(models{i}.cov2-models{j}.cov2)^2/norm(models{i}.cov2)/norm(models{j}.cov2)*10);
			%S(i,j)=exp(-norm(models{i}.cov2-models{j}.cov2,'fro')^2/norm(models{i}.cov2)/norm(models{j}.cov2));
			case 6.6
			S(i,j)=exp(-norm(models{i}.mean-models{j}.mean)^2/eps^2/30)*exp(-norm(models{i}.cov-models{j}.cov,'fro')^2/norm(models{i}.cov+models{j}.cov,'fro')^2);
			case 7
                %if models{i}.dim~=models{j}.dim
                %	S(i,j)=inf;
                %else
                if models{i}.dim==D & models{j}.dim==D
                    S(i,j)=0;
                else
                    %S(i,j)=norm(models{i}.cov2-models{j}.cov2,'fro')^2/((norm(models{i}.cov2,'fro')^2)+(norm(models{j}.cov2,'fro')^2));
                    S(i,j)=norm(models{i}.cov2-models{j}.cov2,'fro')^2/2/max(min(models{i}.dim, size(models{i}.cov2,1) -models{i}.dim),min(models{j}.dim, size(models{j}.cov2,1) -models{j}.dim));%((norm(models{i}.cov2,'fro')^2)+(norm(models{j}.cov2,'fro')^2));
                end
            case 8
                if models{i}.dim~=models{j}.dim
                    S(i,j)=0;
                    %	S(i,j)=1-norm(models{i}.cov2*(models{i}.num-1)-models{j}.cov2*(models{j}.num-1),'fro')^2/((norm(models{i}.cov2*(models{i}.num-1),'fro')^2)+(norm(models{j}.cov2*(models{j}.num-1),'fro')^2));
                    %elseif models{i}.dim==D
                    %	S(i,j)=1;
                else
                    S(i,j)=1-norm(models{i}.cov2*(models{i}.num-1)-models{j}.cov2*(models{j}.num-1),'fro')^2/((norm(models{i}.cov2*(models{i}.num-1),'fro')^2)+(norm(models{j}.cov2*(models{j}.num-1),'fro')^2));
                end
            case 9
                S(i,j)=norm(models{i}.cov2-models{j}.cov2,'fro')^2/((norm(models{i}.cov2,'fro')^2)+(norm(models{j}.cov2,'fro')^2));
            case 10
                S(i,j)=norm(models{i}.cov2*(models{i}.num-1)-models{j}.cov2*(models{j}.num-1),'fro')^2/((norm(models{i}.cov2*(models{i}.num-1),'fro')^2)+(norm(models{j}.cov2*(models{j}.num-1),'fro')^2));
            case 11
                
                S(i,j)=norm(models{i}.cov2-models{j}.cov2,'fro')^2/2/min(models{i}.dim, 1+size(models{i}.cov2,1) -models{i}.dim);
            case 12
                S(i,j)=1-norm(models{i}.cov2*models{j}.cov2,'fro')^2/max((norm(models{i}.cov2,'fro')^2),(norm(models{j}.cov2,'fro')^2));
            case 13
                S(i,j)=det(models{i}.cov)
                
        end
    end
end
%S=S;
if S(1,1)<10^-1
    for i=1:size(S,1)
        S(i,setdiff(1:size(S,2),neighbor{i}))=Inf;
        [a,b]=sort(real(S(i,:)));
        if sum(1./a)==0
            sigmaa(i)=1;
        elseif 1/a(sigma)>0
            sigmaa(i)=1/a(sigma);
        else
            sigmaa(i)=min(1./a(find(1./a>0)));
        end
        
    end
    
    
    sigmaa(find(sigmaa==Inf))=1;
    %S=S.*S1;
    S=real(S);
    
    S1=S;
    %SS=S;
    %for i=1:size(SS,1)
    %SS(i,i)=inf;
    %end
    %SS=reshape(SS,1,size(SS,1)^2);
    %SS=SS(find(SS~=inf));
    %sigmaa=repmat(1/median(SS),1,size(S,1));
    
    %for i=1:length(sigma)
    %S2{i}=exp(-S1.^2/sigma(i)^2);
    %S2{i}=(S2{i}+S2{i}')/2;
    %end
    S2{1}=exp(-repmat((sigmaa),size(S,1),1).*S1.^2.*repmat((sigmaa)',1,size(S,1)));
    S2{1}=(S2{1}+S2{1}')/2;
    clear S1
    S1=S2;
else
    for i=1:size(S,1)
        S(i,setdiff(1:size(S,2),neighbor{i}))=0;
    end
    clear S1
    %S=S./S1;
    S=real(S);
    S1{1}=(S+S')/2;
    for i=1:size(S1{1},1)
        S1{1}(i,i)=1;
    end
end
%for i=1:size(S1,1)
%S1(find(S1(i,:)==Inf))=0;
%end
