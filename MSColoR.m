function [w, v, u, A, B] = MSColoR(X, X_ref, Y, B_logi, paras, t)
% 
% Input:
%          X: genetype data, n * d
%          Y: longitudinal imaging data, g * c (n * m)(number of time points)
%          paras: 1*1
%          paras: g*1
% Output:
%           W: d * c * g, slope values related to gene
%           V: d * c * g, intercept values realated to gene
%           A: 1 * c, values of slope related to aging
%           B: 1 * c, values of intercept related to aging
%           b: d * 1 * g, weights for reg part
 
 % data
[g, c] = size(Y);
[n_1kg, d] = size(X_ref);
m = size(Y{1,1}, 2);

% pre-calculate
S_XX = X_ref' * X_ref / (n_1kg - 1);
e = ones(m, 1);

% initialize
for r = 1:g
    n(r) = size(X{r}, 1);
    w{r} = ones(d, c);
    v{r} = ones(d, c);
    u{r} = ones(d, 1);
end

for j = 1 : c
    for r = 1:g
        for i = 1: n(r) 
            Y_b(i, :) = Y{r, j}(i, :) - X{r}(i, :) * w{r}(:, j) * t(r, :) - X{r}(i, :) * v{r}(:, j); 
        end
        YY1 = mean(Y_b);
        beta(r, :) = polyfit(t(r, :), YY1, 1);
    end
    beta = mean(beta);
    for r = 1:g
        Y_gene{r, j} = Y{r, j} - beta(2) - beta(1) * t(r, :); % r*c(n*m)
    end
end

p = [u', w']; 
q = [u', v']; 

% iteration
tol = 1e-5;

for iter = 1:100

    w_old = w;
    v_old = v;

    for r = 1:g
        XX = X{r}' * X{r};
        XeX = m * XX;
        XttX = sum(t(r, :) .* t(r, :)) * XX;

        % solve w v 
        pl =[p{r, 1}, p{r, 2}];  
        ql = [q{r, 1}, q{r, 2}];
        pi = max(abs(pl), [], 2);  %L∞,1
        Dp = diag(0.5 ./ (pi + eps));
        qi = max(abs(ql), [], 2);
        Dq = diag(0.5 ./ (qi + eps));

        for j = 1 : c
            % update w
            XY = zeros(d, 1);
            for i = 1:n(r)
                XY = XY + (Y_gene{r, j}(i, :) - X{r}(i, :) * v{r}(:, j)) * t(r, :)' * X{r}(i, :)';
            end
            w{r}(:, j) = (XttX + paras(r) * Dp) \ XY;

            % update v
            XY = zeros(d, 1);
            for i = 1:n(r)
                XY = XY + (Y_gene{r, j}(i, :) - X{r}(i, :) * w{r}(:, j) * t(r, :)) * e * X{r}(i, :)';
            end
            v{r}(:, j) = (XeX + paras(r) *Dq)  \ XY;
        end

        % update b
        u{r} = (paras(r) * (Dp +Dq) + S_XX) \ B_logi;
    end
    p = [u', w'];
    q = [u', v'];
    
    % salve beta
	for j = 1 : c
        for r = 1:g
            for i = 1: n(r) 
                Y_b(i, :) = Y{r, j}(i,:) - X{r}(i, :) * w{r}(:, j) * t(r, :) - X{r}(i, :) * v{r}(:, j);     
            end
            YY1 = mean(Y_b);
            beta(r, :) = polyfit(t(r, :), YY1, 1); 
        end
        Beta = mean(beta);
        A(j) = Beta(1);
        B(j) = Beta(2);
        for r = 1:g
            Y_gene{r, j} = Y{r, j} - Beta(2) - Beta(1) * t(r, :); % r*c (n*m)
        end
	end
    
    
    
    % stopping condition
    for r = 1:g
        tol_w(:, r) = max(abs(w{r} - w_old{r}), [], 2);
        tol_v(:, r) = max(abs(v{r} - v_old{r}), [], 2);
    end
    
    if iter == 1
        continue
    end

	if any(tol_w < tol, 'all') || any(tol_v < tol, 'all')
        break
	end
end  




    
 

end



