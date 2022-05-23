function Value=EuroPut(S0,K,r,sigma,T,Smax,n,m)
h=Smax/m;
k=T/n;
% create boundary conditions
x=0:h:Smax;
wxT(1,:)= max(K-x,0);

t=0:k:T;
w0t(:,1)= K*exp(-r*(T-t));
% set up a 
i=1:m-1;
a=(k/2)*(r*i-sigma^2*i.^2);
% set up b
b=1+k*(sigma^2*i.^2+r);
% set up c
c=(k/2)*(-r*i-sigma^2*i.^2);


% set up matrix A, Dimension:(m-1)*(m-1)
A=zeros(m-1);
for i=1:m-1
    A(i,i)=b(i);
end
for i=1:m-2
    A(i+1,i)=a(i+1);
    A(i,i+1)=c(i);
end

% set up boundary conditions
S=zeros(n+1,m+1);
S(n+1,:)=wxT;
S(:,1)=w0t;
S(1,m+1)=0;
% set up the estimates
for j=n:-1:1
    b=[S(j+1,2), S(j+1,3:m-1), S(j+1,m)]-[a(1)*S(j,1), zeros(1,m-3), c(m-1)*S(j,m+1)];
    S(j,2:m)=SolSystemCrout(A,b);
end
% if the initial asset price does not lie on the grid, we must interpolate
% between the two neighboring points.
Value=interp1(x,S(1,:),S0);