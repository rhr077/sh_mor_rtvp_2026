%%% School of Computing and Mathematical Sciences
%%% University of Greenwich, London, United Kingdom
%%% Sheikh Hassan
%%% Email: s.hassan@gre.ac.uk; rhr171@hotmail.com
%%% Modified for Large Scale Block-SOAR from General SOAR

% [2] Hassan, S., Stoyanov, S., Rajaguru, P. and Bailey, C., 2026. Reduced-order
%     modelling for thermal–mechanical analysis of power electronic modules. Finite
%     Elements in Analysis and Design, 253, p.104488.

%%Start_Functions%%

function [Q] = SOAR(M, D, K, N, n, L)
% function [Q,P,T] = BSOAR(M, D, K, N, n, q1)
% Q = SOAR(av, bv, N, n, q1)
%   computes an orthonormal basis Q of the second-order Krylov subspace:
%     span{ Q } = G_k(A,B,v_0).
%   using space-efficient SOAR procedure 
%
% Input:
%   N    :   size of matrices A and B.
%   funav:   a fuction handle y = funav(x) for y = A*x,
%            where x and y are length-n vectors.
%   funbv:   function handle for y = B*x. (see funav)
%   n    :   dimension of second-order Krylov subspace.
%   v0   :   starting vector, of size n-by-1
%
% [Q,P,T] = SOAR(av, bv, N, n, q1)
%   also return the coefficient matrices P and T of the corresponding 
%   compact Krylov decomposition.
%
% References:
%   [1] Yangfeng Su, and Zhaojun Bai, SOAR: a second-order Arnoldi 
%       method for the solution of the quadratic eigenvalue problem,
%       SIAM J. Matrix Anal. Appl., 2005, 26(3): 640-659.
%
% Author:
%   Ding Lu, Fudan Univ. 2014/6/30.
%

input = size(L,2);
n = ceil(n/input);

%Q = zeros(N,n+1);
Q = sparse(N,(n*input)+1);

%P = zeros(N,n+1);
P = sparse(N,(n*input)+1);

%T = zeros(n+1,n+1);
T = sparse((n*input)+1,(n*input)+1);

tol = 1.0*10^-12; % breakdown threshold

% initilize..
temporary_L = distributed(L);
temporary_K = distributed(K);
temporary_q1_dis = temporary_K\temporary_L;
temporary_q1 = gather(temporary_q1_dis);
q1 = sparse(temporary_q1/normest(temporary_q1));
Q(:,1:input) = q1;
f = sparse(N,input);
P(:,1:input) = f;

% temporary 
temporary_M = distributed(M);
temporary_D = distributed(D);
% temporary_K = distributed(K);

for j = 1:(n*input)
    temporary_Q = distributed(Q(:,j));
    temporary_P = distributed(P(:,j));
    temporary_r = (-temporary_K\(temporary_D*temporary_Q + temporary_M*temporary_P));
    r = gather(temporary_r); %Mat-Vec Multiplication
    alpha1 = sparse(normest(r));
    % MGS
    % t = sparse(j,1);
    t = sparse(j+input-1,1); %BSOAR
    for jj = 1:j+input-1
        t(jj) = sparse(Q(:,jj)'*r);
        r = sparse(r - t(jj)*Q(:,jj)); 
    end
    % re-orth
    if normest(r) < 0.717*alpha1 
        for jj = 1:j+input-1
            tt = sparse(Q(:,jj)'*r);
            r = sparse(r - tt*Q(:,jj)); 
            t(jj) = sparse(t(jj) + tt);
        end
    end
    
    T(1:j+input-1,j) = t;
    beta = normest(r);
    T(j+input,j) = beta;
    %T = sparse(T);
    
    % check for breakdown
    if beta>tol
        Q(:,j+input) = sparse(r / beta);
        % Q = sparse(Q);
        ej=sparse(j,1);
        ej(j) = 1;
        f = sparse(Q(:,1:j)*(T(2:j+1,1:j)\ej));
    else
        T(j+1,j) = 1;
        Q(:,j+input) = sparse(N,1);
        ej=sparse(j,1);
        ej(j) = 1;
        f = sparse(Q(:,1:j)*(T(2:j+1,1:j)\ej));
        error('Breakdown in SOAR');
    end
    P(:,j+1) = f;
end

Q = Q(:,1:n*input);
% P = P(:,1:n*input);

end
%END

%%End_Functions%%