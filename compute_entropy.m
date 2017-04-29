%Compute economic bit rate
%Dave Zachariah 2017-04-25

clear

txt = ['dataset2/KOR1990_Tech.csv'];

%% Import flow matrix data
F = dlmread( txt, ';'); %flow matrix
%F = [1 0 0; 0.5 0.5 0; 0.4 0 0]
[n,~] = size(F);
disp(n)

%% Extract basic sector
% Exclude sectors with *all*-zero rows
n_tmp = 0;
while 1

    zvec = zeros(1,n);
    idx = [];
    for i = 1:n
        if sum( F(i,:) == zvec ) ~= n
            idx = union(idx,i);
        end
    end

    %Chip off
    F = F(idx,idx);
    [n_tmp,~] = size(F);
    
    if n == n_tmp
        break
    end
    n = n_tmp;

end


%% Contruct probability matrix
P = zeros(n,n);
for i = 1:n
    P(i,:) = F(i,:) / sum( F(i,:) );
end


%% Compute stationary-probability vector
[V,D]   = eig( P ); %eigendecomposition P = V*D*V'
[~,idx] = max( real(diag(D)) ); %find eigenvalue = 1 (maximum)
pi_vec  = V(:,idx) / sum(V(:,idx)); %compute vector of probabilties (sanity check)


%% Compute entropy rate [in bits]
H = 0;

for i = 1:n
   for j = 1:n

       if P(i,j) > 0
        H = H - ( pi_vec(i) * P(i,j) * log2(P(i,j)) );
       end
       
   end
end

disp(txt)

disp('Number of sectors')
disp(n)

disp('Entropy rate for n sector system [bits per step]')
disp(H)

disp('Fraction of maximum entropy rate [%]')
disp( (H / log2(n)) * 100 )
