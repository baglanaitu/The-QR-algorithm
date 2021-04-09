clear; clc

% Input
n=10; 
A = create_matrix(10);

%% Task 1
A_k = A;
for i=1:10
    [Q,R] = QR(A_k);
    A_k = Q'*A_k*Q;
end
D = diag(A_k)';
D = D(end:-1:1);

disp('Task 1')
disp('________________________')
disp('Calculated eigen vector: ')
disp(D');
disp('Exact eigen vector: ')
disp(exact_values(n)');
disp('Only first eigen values are the same. Then with increasing of')
disp('number of iterations it becomes different as it shown in figure 1.')
disp(' ')


figure (1)
plot(D, 'r-*')
hold on
plot(exact_values(n), 'b-*')
legend('Calculated values', 'Exact values')
title('Eigen values')

%% Task 2

exact_eigen_vector = zeros(n,1);
for k = 1:10
   exact_eigen_vector(k) = sin((pi*k)/(n+1));
end
exact_eigen_vector = exact_eigen_vector / norm(exact_eigen_vector);

% Set up with the smallest value of 
A = A - eye(n)*A_k(n);

disp(' ')
disp('Task 2')
disp('_______________________')
disp('Calculated eigen vector:');
disp(calculate_eigen_vector(n,0.0001,A));
disp('Exact eigen vector:');
disp(exact_eigen_vector);
disp('Using inverse iteration and smallest eigenvalue of A from task 1,')
disp('we could obtain the same results with exact eigenvector values.')

figure (2)
subplot 211
plot(calculate_eigen_vector(n,0.0001,A), 'r-*')
legend('calculated value')
title('Inverse iteration')
subplot 212
plot(exact_eigen_vector, 'b-*');
legend('exact value')


%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of triagonal matrix
function [M] = create_matrix(n)
    M = zeros(n,n);
    for i = 1:n
        for j = 1:n
            if i == j
                M(i,j) = 2;
            elseif j == i-1
                M(i,j) = -1;
            elseif j == i+1
                M(i,j) = -1;
            end
        end
    end
    M = (n+1).^2 .* M;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact values (lambda) 
function [values] = exact_values(n)
    values = zeros(1,n);
    for k=1:n
        values(k) = (4*(n+1)^2) * (sin(k* pi / (2*(n+1))))^2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QR algorithm function
function [Q,R] = QR(M)
    % Initialization
    [n,~] = size(M);
    Q = M;
    Q_s = zeros(n,n,n);
    for k=1:n
        e = double(zeros(n-k+1,1));
        % Get lower part of current column (from the diagonal) 
        lower_part = Q(k:end,k);
        e(1,1) = norm(lower_part);
        w = (lower_part + e*sign(lower_part(1,1)));
        v = w./norm(w);
        
        Q_k = diag(ones(1,n));
        Q_k(k:end,k:end) = eye(size(v*v'))-2*(v*v');
        Q_s(:,:,k) = Q_k;
        Q = Q_k*Q;
    end
    R = Q;
    % Recompose the final Q matrix with combination of intermediary matrices
    Q = diag(ones(1,n));
    for i=1:n
       Q = Q*Q_s(:,:,i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v] = calculate_eigen_vector(n, threshold,A)
    v = 1/n .* ones(n,1);
    % ==> norm v = 1

    v_previous = zeros(n,1);
    
    while( norm(v_previous-v) > threshold)
        v_previous = v;
        % Solve Aw=v
        w = linsolve(A,v); 
        v = w/norm(w);
    end
end