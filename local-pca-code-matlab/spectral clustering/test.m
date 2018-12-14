N=3000;
theta=[0:1/N:(N-1)/N].^2*pi*2;
theta=theta.*sign(rand(1,N)-0.5);
X=[cos(theta);sin(theta)];

saveas(plot(X(1,:),X(2,:),'.'),'temp.fig')


Y=repmat(sum(X.^2,1),N,1)+repmat(sum(X.^2,1),N,1)'-2*X'*X;
Y=exp(-Y/0.02);
sumY=sum(Y,2);
Z=diag(sumY.^-1)*Y*diag(sumY.^-1);
sumZ=sum(Z,2);
Z=diag(sumZ.^-0.5)*Z*diag(sumZ.^-0.5);

[U,S]=svd(Z);
diag(S(1:10,1:10))
X1=U(:,1:2)*(S(1:2,1:2).^0.5);
Y1=zeros(size(X1));
U1=Y1;
for i=1:N
Y1(i,:)=X1(i,:)/norm(X1(i,:));
end
for i=1:N
U1(i,:)=U(i,1:2)/norm(U(i,1:2));
end
saveas(plot(X1(:,1),X1(:,2),'.'),'temp1.fig')
saveas(plot(Y1(:,1),Y1(:,2),'.'),'temp2.fig')
saveas(plot(U(:,1),U(:,2),'.'),'temp3.fig')
saveas(plot(U1(:,1),U1(:,2),'.'),'temp4.fig')
