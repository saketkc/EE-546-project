N1=100;N2=100,;N3=100;X=rand(N1,2);X1=rand(N2,2)*0.02+repmat([0.2,0.2],N2,1);X1=rand(N3,2)*0.02+repmat([0.8,0.8],N3,1);M=[ones(1,N1),2*ones(1,N2),3*ones(1,N3)];X=[X;X1;X2];