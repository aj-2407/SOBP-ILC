function [ Q R ] = QR_decompos( X )
%名称：QR分解
%功能：将矩阵X进行QR分解（按Gram_Schmidt法则正交化），返回分解后的矩阵QR
%Date：20160601
%调用实例：A = [1 0 0;1 1 0;1 1 1;1 1 1];[Q R] = QR_decompos(A)

%Begin-------------------------------------
[m n] = size(X); %m行数，n列数
Q = zeros(m,n);
R = zeros(n,n);

%容错控制
if m<n
    error('<<矩阵行数小于列数，无法保证列满秩,不能进行QR分解！>>')
end

%正交化
Q(:,1) = X(:,1);
for i = 2:n
    sum = zeros(m,1);
    for j=1:i-1
        sum = sum + X(:,i)'*Q(:,j)/(Q(:,j)'*Q(:,j))*Q(:,j);
    end
    Q(:,i) = X(:,i) - sum;  
end

%单位化
for i = 1:n
    Q(:,i) = Q(:,i)/norm(Q(:,i));  
end

%求R
R = Q'*X;

%End---------------------------------------
end

