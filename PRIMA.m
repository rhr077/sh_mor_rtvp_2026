%%% School of Computing and Mathematical Sciences
%%% University of Greenwich, London, United Kingdom
%%% Sheikh Hassan
%%% Email: s.hassan@gre.ac.uk; rhr171@hotmail.com
%%% Modified for Large Scale PRIMA from General PRIMA

% [1] Hassan, S., Rajaguru, P., Stoyanov, S., Bailey, C. and Tilford, T., 2024.
%     Coupled thermal-mechanical analysis of power electronic modules with finite
%     element method and parametric model order reduction. Power Electronic Devices
%     and Components, 8, p.100063.
% [2] Hassan, S., Stoyanov, S., Rajaguru, P. and Bailey, C., 2026. Reduced-order
%     modelling for thermal–mechanical analysis of power electronic modules. Finite
%     Elements in Analysis and Design, 253, p.104488.

%% Start Function
function [V]=PRIMA(E,A,B,order)

% PRIMA ORDER REDUCTION
% Reduces a linear time invariant system of the form
%
% E dx/dt =A x+ B u
% y = C x
%
% to an approximate system of the desired order
% using the PRIMA algorithm.
%

% E: n x n
% A: n x n
% B: n x ni, where ni is the number of inputs
% C: n x no, where no is the number of outputs
% r: desired order for the reduced model (integer)
%
%
% V: n x order projection matrix
%
% Process input arguments and options
% Out put V and W are the projection matrices
% Reduced matrices
% Er = V'*E*V
% Ar = V'*A*V
% Br = V'*B
% Cr = C*V

s0=0;

numorth=2;

% Order of original system

n=size(A,1);

% Number of inputs and outputs

% ni=size(R,2);
ni=size(B,2);

% no=size(Q,2);


if order >= n

error('The required order is greater or equal to the order of the original system!');

end

% Number of moments to be matched
k=ceil(order/ni);

% Aq=zeros(k*ni,k*ni);
Aq=sparse(k*ni,k*ni);

% Vq=zeros(n,k*ni);
Vq=sparse(n,k*ni);

% delta=zeros(ni,ni);
% delta=sparse(ni,ni);

% I=eye(n);
% I = speye(n);
temporary_E = distributed(E);
temporary_A = distributed(A);
temporary_B = distributed(B);
temporary_R = (s0*temporary_E + temporary_A)\temporary_B;
R = gather(temporary_R); 

% Generate first block V_0 of projection matrix
% R = gather(R);
[Vq(:,1:ni),~] = qr(R,0);

% Arnoldi iteration

for j=1:k-1

    % % Vq(:,j*ni+1:(j+1)*ni) = -(G+s0*C)\(C*Vq(:,(j-1)*ni+1:j*ni));
    % 
    temp = Vq(:,(j-1)*ni+1:j*ni);
    
    temp = distributed(temp);

    % V = E * Vq(:,j*ni+1:(j+1)*ni);
    % Vq(:,j*ni+1:(j+1)*ni) = A\V;
    % Vq(:,j*ni+1:(j+1)*ni) = A\(E*Vq(:,j*ni+1:(j+1)*ni));
    % temp = A\(E*temp);
	
    temp = temporary_A\(temporary_E*temp); % A\(E*V) for PRIMA
    temp = gather(temp);
    Vq(:,j*ni+1:(j+1)*ni) = temp;

    for xi=1:numorth

        for i=1:j % Modified Gram-Schmidt orthonormalisation

            delta = Vq(:,(j-i)*ni+1:(j-i+1)*ni)' * Vq(:,j*ni+1:(j+1)*ni);

            Vq(:,j*ni+1:(j+1)*ni) = Vq(:,j*ni+1:(j+1)*ni) - Vq(:,(j-i)*ni+1:(j-i+1)*ni) * delta;

        end

    end
    [Vq(:,j*ni+1:(j+1)*ni),Aq(j*ni+1:(j+1)*ni,(j-1)*ni+1:j*ni)] = qr(Vq(:,j*ni+1:(j+1)*ni),0);
end

% Trim Vq in order to achieve the desired order

Vq=Vq(:,1:order);

% Assign the projection matrix V

V=Vq;

% return

end

%%End_Functions%%