%%% School of Computing and Mathematical Sciences
%%% University of Greenwich, London, United Kingdom
%%% Sheikh Hassan
%%% Email: s.hassan@gre.ac.uk; rhr171@hotmail.com
%%% Modified for Large Scale Block-Arnoldi from General Arnoldi

% [1] Hassan, S., Rajaguru, P., Stoyanov, S., Bailey, C. and Tilford, T., 2024.
%     Coupled thermal-mechanical analysis of power electronic modules with finite
%     element method and parametric model order reduction. Power Electronic Devices
%     and Components, 8, p.100063.
% [2] Hassan, S., Stoyanov, S., Rajaguru, P. and Bailey, C., 2026. Reduced-order
%     modelling for thermal–mechanical analysis of power electronic modules. Finite
%     Elements in Analysis and Design, 253, p.104488.

%%Start_Functions%%

function [V]=Arnoldi(E,A,B,r)

% This method computes the projection matrices V and W
% by Arnoldi algorthm (AA) MOR technique
%
% E dx/dt =A x+ B u
% y = C x
%
% to an approximate system of the desired order
% using the Arnoldi algorithm.
%
% E: n x n
% A: n x n
% B: n x ni, where ni is the number of inputs
% C: n x no, where no is the number of outputs
% r: desired order for the reduced model (integer)
%

% V,W: n x order projection matrix
temporary_E = distributed(E);
temporary_A = distributed(A);
temporary_B = distributed(B); 
temporary_C = temporary_A\temporary_B;
C = gather(temporary_C);
input = size(C,2);
k = ceil(r/input);
N = size(A,1);
V = sparse(N,k*input);

numort = 2;

% C = gather(C);
[temp,~] = qr(C,0);
% V(:,1)=(A\B)/sqrt((A\B)'*(A\B));
% V(:,1)=(C)/sqrt((C)'*(C));
V(:,1:input) = temp;
clear temp;
% V = gather(V); 

% W=(A'\C')/sqrt((A'\C')'*(A'\C'));

for i=1:k-1
    
    % V = distributed(V);
    temp = V(:,((i-1)*input)+1:i*input);
    temp = distributed(temp);
    temp = temporary_E\(temporary_A*temp); % E\A*V; for Arnoldi 
    temp = gather(temp);
    V(:,(i*input)+1:(i+1)*input) = temp;
    clear temp;
    % V(:,i)=E\(A*V(:,i-1));
    % V = gather(V);
    
    % W(:,i)=A'\W(:,i-1);
    
    for xi=1:numort

        for j=1:i
    
            h1 = V(:,((i-j)*input)+1:(i-j+1)*input)' * V(:,(i*input)+1:(i+1)*input);
        
            % h2=W(:,i)'*W(:,j);
        
            V(:,(i*input)+1:(i+1)*input) = V(:,(i*input)+1:(i+1)*input) - V(:,((i-j)*input)+1:(i-j+1)*input) * h1;
        
            % W(:,i)=W(:,i)-h2*W(:,j);

        end

    end

    temp = V(:,(i*input)+1:(i+1)*input);
    [temp,~] = qr(temp,0);
    V(:,(i*input)+1:(i+1)*input) = temp;
    % clear temp;
    % V(:,i)=V(:,i)/sqrt(V(:,i)'*V(:,i));
    
    % W(:,i)=W(:,i)/sqrt(W(:,i)'*W(:,i));

end

V = V(:,1:r);

% return

end  