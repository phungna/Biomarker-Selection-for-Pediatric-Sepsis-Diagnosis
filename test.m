clear all

deg_genes = readtable("ReGeneL\dulieudaxuly.csv")
deg_genes(:,{'Var1', 'SID'}) = []

%split data
X = deg_genes(:, 1:108)
y = deg_genes(: , 109)

%normalize data
x  = X.Variables
for i = 1:108
    x(:, i) = normalized(x(:, i))
end
x
Y = y.Variables

%train test spit
cv = cvpartition(height(deg_genes), 'HoldOut', 0.50)
x_train = x(cv.training, :);
y_train = Y(cv.training, :);
x_test = x(cv.test, :);
y_test = Y(cv.test, :);
