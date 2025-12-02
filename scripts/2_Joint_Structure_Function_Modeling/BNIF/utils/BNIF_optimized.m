function [f, P] = BNIF(D, R, varargin)

% BNIF - Optimized Brain Network Information Flow mapping
% Uses IBM CPLEX solver if available, otherwise falls back to linprog
% Inputs:
%   D - NxN adjacency matrix
%   R - NxM matrix of functional activation modes
%   varargin - optional parameters: rho, r, MaxIter, SolnTH
% Outputs:
%   f - NxNxM matrix of functional flows
%   P - NxN matrix of connectivity adjustments
%
%% Input Validation
ip = inputParser;
ip.CaseSensitive = true;
addRequired(ip, 'D', @(x) isnumeric(x) && all(x(:) >= 0) && issymmetric(x));
addRequired(ip, 'R', @(x) isnumeric(x) && all(x(:) >= 0));
addOptional(ip, 'rho', 0.01, @(x) isnumeric(x) && isscalar(x) && (x >= 0));
addOptional(ip, 'r', 0, @(x) isnumeric(x) && isscalar(x) && (x >= 0));
addParameter(ip, 'MaxIter', int32(1e5), @(x) isscalar(x) && (x > 0) && isinteger(x));
addParameter(ip, 'SolnTH', 1e-6, @(x) isnumeric(x) && isscalar(x) && (x >= 0));
parse(ip, D, R, varargin{:});
if size(D,1) ~= size(R,1)
    error('D and R must have the same number of rows');
end

MODEL.R = R;
[MODEL.N, MODEL.M] = size(R);
MODEL.L = (MODEL.N^2 - MODEL.N) / 2;
MODEL.D = D(triu(true(MODEL.N), 1));
MODEL.rho = ip.Results.rho;
if ip.Results.r == 0
    R_mode = sum(R, 2);
    D_node = sum(D, 2);
    r = R_mode ./ D_node;
    MODEL.r = max(r(D_node > 0));
else
    MODEL.r = ip.Results.r;
end
M = MODEL.M; L = MODEL.L; N = MODEL.N;

%% Phase 0: Linear Program to Estimate P
LP0.lb = sparse(M*L + L, 1);
LP0.ub = Inf * ones(M*L + L, 1);
LP0.b = sparse(L + N*M + L*M, 1);
LP0.A = sparse(L + N*M + L*M, M*L + L);

rows = repmat((1:L)', M, 1);
cols = reshape(bsxfun(@plus, (0:M-1)*L, (1:L)'), [], 1);
vals = ones(L*M, 1);
A_block = sparse(rows, cols, vals, L, L*M);
LP0.A(1:L, 1:L*M) = A_block;
LP0.A(1:L, L*M+1:L*M+L) = -MODEL.r * speye(L);
LP0.b(1:L) = MODEL.r * MODEL.D;

[row_i, row_j] = find(triu(true(N),1));
link_map = containers.Map;
for l = 1:L
    link_map(sprintf('%d_%d', row_i(l), row_j(l))) = l;
end

row_idx = cell(M,1); col_idx = cell(M,1); val_idx = cell(M,1); b_vals = cell(M,1);
parfor m = 1:M
    r_idx = []; c_idx = []; v_idx = []; b_m = zeros(N,1);
    for n = 1:N
        idx = L + (m-1)*N + n;
        for j = 1:N
            if j == n, continue; end
            i1 = min(n, j); i2 = max(n, j);
            key = sprintf('%d_%d', i1, i2);
            if isKey(link_map, key)
                l_idx = link_map(key);
                r_idx(end+1) = idx;
                c_idx(end+1) = (m-1)*L + l_idx;
                v_idx(end+1) = -1;
            end
        end
        b_m(n) = -MODEL.R(n, m);
    end
    row_idx{m} = r_idx; col_idx{m} = c_idx; val_idx{m} = v_idx; b_vals{m} = b_m;
end
LP0.A = LP0.A + sparse(vertcat(row_idx{:}), vertcat(col_idx{:}), vertcat(val_idx{:}), L+N*M+L*M, L*M+L);
LP0.b(L+1:L+N*M) = vertcat(b_vals{:});

row_idx = (1:L*M)';
col_idx = (1:L*M)';
val_idx = ones(L*M, 1);
LP0.A(L+N*M+1:end, 1:L*M) = sparse(row_idx, col_idx, val_idx, L*M, L*M);

b_idx = L+N*M;
b_parts = zeros(L, M);
parfor m = 1:M
    l = 1; b_vec = zeros(L,1);
    for i = 2:N
        for j = 1:i-1
            if any([MODEL.R(i,m), MODEL.R(j,m)] == 0)
                b_val = 0;
            else
                b_val = max(MODEL.R(i,m), MODEL.R(j,m));
            end
            b_vec(l) = b_val;
            l = l + 1;
        end
    end
    b_parts(:, m) = b_vec;
end
for m = 1:M
    LP0.b(b_idx + (m-1)*L + (1:L)) = b_parts(:, m);
end

minD = min(MODEL.D(MODEL.D > 0));
invD = 1 ./ (MODEL.D + min(1e-6, minD)*(MODEL.D == 0));
LP0.c = sparse([repmat(invD, M, 1); MODEL.rho * (1 + invD)]);

useCPLEX = exist('cplexlp', 'file') == 2 || exist('cplexlp', 'file') == 6;
if useCPLEX
    disp('Using CPLEX for Phase 0');
    try
        cplexopts = cplexoptimset('cplex');
        cplexopts.threads = feature('numCores');
	disp(cplexopts.threads);
    catch
        warning('Could not configure CPLEX thread count, using default');
        cplexopts = [];
    end
    [opt_x, ~, EXITFLAG] = cplexlp(LP0.c, LP0.A, LP0.b, [], [], LP0.lb, LP0.ub);
else
    disp('Using linprog for Phase 0');
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'MaxIter', ip.Results.MaxIter, 'Display', 'none');
    [opt_x, ~, EXITFLAG] = linprog(LP0.c, LP0.A, LP0.b, [], [], LP0.lb, LP0.ub, options);
end
if EXITFLAG <= 0
    error('Phase 0 optimization failed. Exitflag: %d', EXITFLAG);
end
opt_x = opt_x .* (opt_x > ip.Results.SolnTH);

Pl = zeros(N);
Pl(triu(true(N),1)) = opt_x(M*L+1:end);
Pl = Pl + Pl';

%% Phase 1: Minimize max f / (D + P)
DP = MODEL.D + Pl(triu(true(N),1));
LP1_lb = LP0.lb(1:M*L);
LP1_ub = LP0.ub(1:M*L);
LP1_lb(end+1) = 0;
LP1_ub(end+1) = 2 * MODEL.r;
LP1_b = LP0.b;
LP1_b(1:L) = 0;
LP1_A = [LP0.A(:,1:M*L), sparse(size(LP0.A,1),1)];
LP1_A(1:L, M*L+1) = -DP;
LP1_c = sparse(M*L+1,1); LP1_c(end) = 1;

if useCPLEX
    disp('Using CPLEX for Phase 1');
    [opt_x, ~, EXITFLAG] = cplexlp(LP1_c, LP1_A, LP1_b, [], [], LP1_lb, LP1_ub);
else
    disp('Using linprog for Phase 1');
    [opt_x, ~, EXITFLAG] = linprog(LP1_c, LP1_A, LP1_b, [], [], LP1_lb, LP1_ub, options);
end
if EXITFLAG <= 0
    error('Phase 1 optimization failed. Exitflag: %d', EXITFLAG);
end
opt_x = opt_x .* (opt_x > ip.Results.SolnTH);

f = zeros(N, N, M);
parfor m = 1:M
    tmp = zeros(N);
    tmp(triu(true(N),1)) = opt_x((m-1)*L+1:m*L);
    f(:,:,m) = tmp + tmp';
end

P = Pl;
end
