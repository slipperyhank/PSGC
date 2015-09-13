function [L,Gamma]=PSGC_noIC(dN,R)

% Variables for error recording, not finished yet

flag=0;

global delta


% convergence criteria
epsilon=0.00001;

% number of segments
K=length(dN);

%number of parameters
S=length(R(:)) / K;

%Initialize output
Gamma=zeros(S,1);
L=0;





% Initial Conditions
X=zeros(S,1);
%size(X)
index=1:S;

% Bad measures parameters which go below the minimum (-10)
Bad=X;


%initialize Jacobian (J) and Derivative (F) and intensity (lambda)
J=zeros(S);
F=zeros(S,1);
lambda=zeros(K,1);

while 1==1
    for k=1:K
        % Calculate lambda
        size(X);
        size(R(k,:));
        lambda(k)=exp(sum(X'.*R(k,:))); 
        % If probability of a shift grows too large, model fails
        if lambda(k)*delta > 1
            flag=1;
            break
        end
    end
    if flag==1
        break
    end
    
    % derivative (F) and jacobian (J) of parameters
    for k=1:length(X)
        F(k) = sum((dN.*R(:,k) + R(:,k).*(delta*lambda.*(dN-1))./(1-delta.*lambda)));
    end
    for k=1:length(X)
        for j=1:length(X)
            J(k,j) = sum((R(:,k).*R(:,j).*(dN-1).*delta.*lambda./(1-lambda.*delta).^2));
        end
    end
    % Calculate change in parameters dX. Magnitude is reduced to 20% to
    % prevent parameters from over shooting and jumping outside the support 
    % of the model
    F=-F;
    dX = J\F;
    X = X + dX/5;
    
    % If any of the parameters are below -10, go to recursion step
    if min(X) < -10
        Bad=X<-10;
        flag=2;
        break;
    end
    % when convergence is reached, exit loop
    if max(abs(dX)) < epsilon 
        break
    end
end

% flag=0: No errors, calculate output
if flag==0
    Gamma=X;
    for k=1:K
        lambda(k)=exp(sum(X'.*R(k,:))); 
    end
    L = sum((dN.*log(lambda.*delta) + (1-dN).*log(1-lambda.*delta)));
end

% flag=1: error, must be looked at
if flag==1
    L=0;
    Gamma=X.*0;
end

% flag=2: parameters below -10, remove from optimization and called
% function again
if flag==2
    Gamma(index(Bad==1))=-10;
    index(Bad==1)='';
    ROld='';
    %size(X(Bad==0))
    [L,temp]=PSGC2(dN,R,Bad,ROld,X(Bad==0));
    Gamma(index)=temp;
end
