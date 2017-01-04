function [partition, conductance] = Low_Conductance(A, D, Index, lambda)

%A = A*diag(D)^(-1)*A;
OD = D;

D = D + 0*lambda*ones(1,length(D));

TriA = triu(A);


%MIP GUROBI
u = 10000;
l = 1e-5;
n = length(A);

% construct linear constarint
% cosntraint layout
% [z z*x1 z*x2 ... z*xn x1 x2 ... xn x1*x2 ...]
% linearization
n_var = 1 + 2*n + sum(sum(TriA>0.0));
line_count = 1;

for i = 1 : n
    %z -zx1 + u*x1 <= u
    W(line_count,1) = 1;
    W(line_count,i+1) = -1;
    W(line_count,i+1+n) = u;
    b(line_count) = u;
    line_count = line_count + 1;
    
    % -zx1 + l*x1 <=0
    W(line_count,i+1) = -1;
    W(line_count,i+1+n) = l;
    b(line_count) = 0;
    line_count = line_count + 1;
    
    
%     % -z + zx1 - l*x1 <= -l
%     W(line_count,1) = -1;
%     W(line_count,i+1) = 1;
%     W(line_count,i+1+n) = -l;
%     b(line_count) = -l;
%     line_count = line_count + 1;
%     
%     % zx1 - u*x1 <= 0
%     W(line_count,1+i) = 1;
%     W(line_count,i+1+n) = -u;
%     b(line_count) = 0;
%     line_count = line_count + 1;
end

count = 1;
w = [];
for i = 1 : length(A)
    for j = i+1 : length(A)
        if(A(i,j)>0.0)
            w = [w A(i,j)];
            % -xi + xij <= 0
            W(line_count,i+1+n) = -1;
            W(line_count,count+2*n+1) = 1;
            b(line_count) = 0;
            line_count = line_count + 1;
            
            % -xj + xij <= 0
            W(line_count,j+1+n) = -1;
            W(line_count,count+2*n+1) = 1;
            b(line_count) = 0;
            line_count = line_count + 1;
            
            % xi + xj - xij <= 1
            W(line_count,i+1+n) = 1;
            W(line_count,j+1+n) = 1;
            W(line_count,count+2*n+1) = -1;
            b(line_count) = 1;
            line_count = line_count + 1;
            
            % -xij <=0
            W(line_count,count+2*n+1) = -1;
            b(line_count) = 0;
            line_count = line_count + 1;
            
            count = count + 1;
        end
    end
end

W(line_count, 2:1+n) = D;
W(line_count, 1+2*n+1:end) = -2*w;
b(line_count) = 0;
line_count = line_count+1;

W(line_count, 1+n+1:1+2*n) = 1;
b(line_count) = n;
line_count = line_count+1;

W(line_count, 1+n+1:1+2*n) = -1;
b(line_count) = -3;
line_count = line_count+1;

% W(line_count, 1+n+1) = 1;
% b(line_count) = 1;
% line_count = line_count+1;
% 
% W(line_count, 1+n+1) = -1;
% b(line_count) = -1;



n_var = length(W(1,:));
c = zeros(n_var,1);
c(1) = 1;
objtype = -1; % 1 - minimize, -1 - maximize
W = sparse(W);
b = b';
lb = [];
ub = [];
contypes = '<';
vtypes = [];
for i = 1 : n_var
    vtypes = [vtypes 'C'];
end
for i = 1+n+1 : 1+2*n
    vtypes(i) = 'B';
end

clear opts
opts.IterationLimit = 1000000;
opts.FeasibilityTol = 1e-6;
opts.IntFeasTol = 1e-5;
opts.OptimalityTol = 1e-6;
opts.Method = 1; % 0 - primal, 1 - dual
opts.Presolve = 2; % -1 - auto, 0 - no, 1 - conserv, 2 - aggressive
opts.Display = 1;
opts.LogFile = 'test.log';
opts.WriteToFile = 'test.mps';
[x,val,exitflag,output] = gurobi_mex(c,objtype,W,b,contypes,lb,ub,vtypes,opts);
% Gurobi gives no Lagrange multipliers or reduced costs for MIPs, without calling modelfix

%disp('Solution:');disp(x')
disp('Optimal obj value:');disp(val)
disp('Exit flag:');disp(exitflag)
%disp('Optimization info:');disp(output)

if(exitflag~=2)
    partition = [];
    conductance = 0;
else
    ind = find(x(1+n+1:1+2*n)>0);
    partition = Index(ind);
    conductance = sum(sum(A(ind,ind)))/sum(OD(ind));
end



















