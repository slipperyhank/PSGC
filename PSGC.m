function [L,Gamma]=PSGC(dN,R,x0)


global delta


% convergence criteria
epsilon=0.00001;

% number of segments
K=length(dN);

%number of parameters
S=length(R(:)) / K;
index=1:S;
lambda=zeros(K,1);
flag=0;

%Initial Conditions
X=x0;

if sum(X<=-10) == 0
    
    Bad=X;
    index=1:S;
    %initialize Jacobian (J) and Derivative (F)
    J=zeros(S);
    F=zeros(S,1);

    while 1==1
        for k=1:K
            lambda(k)=exp(sum(X'.*R(k,:))); 
            if lambda(k)*delta > 1
                flag=1;
                break
            end
        end
        if flag==1
            break
        end
    
        for k=1:length(X)
            F(k) = sum((dN.*R(:,k) + R(:,k).*(delta*lambda.*(dN-1))./(1-delta.*lambda)));
        end
        for k=1:length(X)
            for j=1:length(X)
                J(k,j) = sum((R(:,k).*R(:,j).*(dN-1).*delta.*lambda./(1-lambda.*delta).^2));
            end
        end
    
        F=-F;
        dX = J\F;
    
        X = X + dX/5;
    
    
        if min(X) <= -10
            Bad=X<=-10;
            flag=2;
            break;
        end
        if max(abs(dX)) < epsilon 
            break
        end
    end
else
    Bad=X<=-10;
    flag=2;
end
if flag==0
    Gamma=X;
    for k=1:K
        lambda(k)=exp(sum(X'.*R(k,:))); 
    end
    L = sum((dN.*log(lambda.*delta) + (1-dN).*log(1-lambda.*delta)));
end

if flag==1
    L=0;
    Gamma=X*0;
end

if flag==2
    Gamma(index(Bad==1))=-10;
    index(Bad==1)='';
    ROld='';
    [L,temp]=PSGC2(dN,R,Bad,ROld,X(Bad==0));
    Gamma(index)=temp;
end
