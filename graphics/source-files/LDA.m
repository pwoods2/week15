function [a_n,b_n,resub] = LDA( x0,x1 )

% compute the LDA parameters for x0,x1, where the rows of x0 are different
% samples and the column of x0 are the different features(genes).

size_x0 = size(x0,1);
size_x1 = size(x1,1);
p_x1  = size_x1/(size_x0+size_x1);  % prior of class 1
u0 = mean(x0);
u1 = mean(x1);
u_diff = (u1-u0)';
u_sum = (u0+u1)';
cov0 = cov(x0);
cov1 = cov(x1);

inv_cov0 = inv(cov0);
inv_cov1 = inv(cov1);

cov_pooled = ((size_x0-1)*cov0+(size_x1-1)*cov1)/(size_x0+size_x1-2);

inv_cov_pooled = inv(cov_pooled);


a_n = inv_cov_pooled*u_diff;
b_n = -0.5*a_n'*u_sum + log(p_x1/(1-p_x1));


resub = 0;

for i = 1:size_x0
    if((a_n'*x0(i,:)'+b_n)>0)
        resub = resub+1;
    end
end

for i = 1:size_x1
    if((a_n'*x1(i,:)'+b_n)<0)
        resub = resub+1;
    end
end

resub = resub/(size_x0+size_x1);


