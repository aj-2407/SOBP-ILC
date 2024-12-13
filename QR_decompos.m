function [ Q R ] = QR_decompos( X )
%���ƣ�QR�ֽ�
%���ܣ�������X����QR�ֽ⣨��Gram_Schmidt�����������������طֽ��ľ���QR
%Date��20160601
%����ʵ����A = [1 0 0;1 1 0;1 1 1;1 1 1];[Q R] = QR_decompos(A)

%Begin-------------------------------------
[m n] = size(X); %m������n����
Q = zeros(m,n);
R = zeros(n,n);

%�ݴ����
if m<n
    error('<<��������С���������޷���֤������,���ܽ���QR�ֽ⣡>>')
end

%������
Q(:,1) = X(:,1);
for i = 2:n
    sum = zeros(m,1);
    for j=1:i-1
        sum = sum + X(:,i)'*Q(:,j)/(Q(:,j)'*Q(:,j))*Q(:,j);
    end
    Q(:,i) = X(:,i) - sum;  
end

%��λ��
for i = 1:n
    Q(:,i) = Q(:,i)/norm(Q(:,i));  
end

%��R
R = Q'*X;

%End---------------------------------------
end

