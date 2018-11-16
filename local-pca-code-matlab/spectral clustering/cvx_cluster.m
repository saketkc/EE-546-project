function X=cvx_cluster(A);
[N,K]=size(A);
cvx_begin quiet
    variable X(N,K) ;
	minimize (trace(X'*A))
	subject to 
	X>=0;
	X*ones(K,1)==1;
	ones(1,N)*X>= floor(0.1*N)+1;
cvx_end