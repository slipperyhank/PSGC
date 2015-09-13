function [L,Gamma]=PSGC2(dN,R,Bad2,ROld,x0)

global delta

R1=R(:,Bad2==0);
R2=R(:,Bad2==1);
R2=[ROld,R2];


K=length(dN);
S1=length(R1(:))/K;
S2=length(R2(:))/K;
index=1:S1;
% convergence criteria
epsilon=0.00001;

lambda=zeros(K,1);
flag=0;
%Initial Conditions

X=x0;

Y=zeros(S2,1)-10;


Bad=X;

J=zeros(S1);
F=zeros(S1,1);

while 1==1
    for k=1:K
        lambda(k)=exp(sum(X'.*R1(k,:)) + sum(Y'.*R2(k,:))); 
        if lambda(k)*delta > 1
            flag=1;
            break
        end
    end
    if flag==1
        break
    end
    for k=1:S1
        F(k) = sum((dN.*R1(:,k) + R1(:,k).*(delta*lambda.*(dN-1))./(1-delta.*lambda)));
    end
    for k=1:length(X)
        for j=1:length(X)
            J(k,j) = sum((R1(:,k).*R1(:,j).*(dN-1).*delta.*lambda./(1-lambda.*delta).^2));
        end
    end
    F=-F;
    dX = J\F;
    X = X + dX/5;
    
    if min(X) < -10
        Bad=X<-10;
        flag=2;
        break;
    end
    if max(abs(dX)) < epsilon
        break
    end

end

if flag==0
    Gamma=X;
    for k=1:K
        lambda(k)=exp(sum(X'.*R1(k,:)) + sum(Y'.*R2(k,:)));
    end
    L = sum((dN.*log(lambda.*delta) + (1-dN).*log(1-lambda.*delta)));
end

if flag==1
    L=0;
    Gamma=X.*0;
end

if flag==2
    Gamma(index(Bad==1))=-10;
    index(Bad==1)='';
    [L,temp]=PSGC2(dN,R1,Bad,R2,X(Bad==0));
    Gamma(index)=temp;
end
