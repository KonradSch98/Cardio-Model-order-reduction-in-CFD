%% this script is a simple implementation of the DMD algorithm for large 
%% data sets for a more memory efficient calculation

cd 'D:\Eigene Dokumente\Uni\BA\BA_cardioMOR\Testdaten\Dynamisch\Unsteady_1_cycles'
dt = 0.1;
t = 10;

tic

%preperations
folder_content = dir();
X = zeros(5136636,10);

%reads every file seperately. should be adapted to the compact load of the
%other scripts
k=1;
for i = 5:14
    %open and read file
    fileid = folder_content(i).name;
    data = readmatrix(fileid);
    snap = reshape( data(:,7:9) , [] , 1 );
    %store data in snapshot matrices
    X(:,k) = snap;
    k=k+1;
end
clear data


%calculate H matrix
H = zeros(9,9);
for i =1:9
    for j=1:9
        if i<=j
            H(i,j) = dot( X(:,i) , X(:,j) );
            if i ~= j
                H(j,i) = H(i,j);
            end
        end
    end
end

%eigenddecomposition of H
[V,S2] = eig( H );
s = diag(S2);
%sort in descending order like done in SVD
[v , in] = sort(s, 'descend');
V = V(:,in);
S2 = diag( s(in) );
%define additional H matrices
H_p = H(:,2:9);
H_pp = zeros(9,1);
for i = 1:9
    H_pp(i) = dot( X(:,10), X(:,i) );
end

%compute approximation A and its eigenvalues
A_hat = S2^(-0.5) * V' * [H_p, H_pp] * V *  S2^(-0.5);
[W,Lamb1] = eig(A_hat);

%rescale eigenvalues to original matrix
omega = diag(Lamb1)/dt;

%compute scaling value
d = eye(9)/(W' * W) * W' * S2^(-0.5) * V' * H(:,1);
D = diag(d);


%compute the eigenvectors /  dmd_modes 'Phi'
Phi = zeros(5136636,9);
T= V * S2^(-0.5) * W *D;
for i = 1:9
    phi = zeros(5136636,1);
    phi = phi +  X(:,1) * T(i,1);
    Phi(:,i) = phi;
end
