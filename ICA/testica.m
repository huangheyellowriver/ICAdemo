% Parameters of the toy example.
n = 200;    % number of samples (quick)
%n = 1000;  % number of samples (precise)
sd1 = 0.15; % standard deviation of source 1
sd2 = 0.8;  % standard deviation of source 2
angle = 30; % rotation angle


% Generate the bimodal data set.
randn('state',1);
s1 = (sd1 * randn(1,n) + sign(randn(1,n)));
s2 = sd2 * randn(1,n);
S = [s1; s2];

% Generate the mixing matrix to be a rotation.
theta = angle / 360 * 2*pi;
A = [cos(theta) sin(theta); ...
    -sin(theta) cos(theta)];

% Linearly mix the data.
X = A * S;

% Subtract off the mean of each dimension.
X = X - repmat(mean(X,2),1,n);

% Calculate the whitening filter.
[E, D] = eig(cov(X'));

% Whiten the data
X_w = sqrtm(pinv(D))*E'*X;

% Create an array of angles between 0 and 180 deg.
angles = 0:2:180;
thetas = angles / 360 * 2*pi;

% Calculate the multi-information for all angles.
for i=1:length(thetas)
    % Generate a rotation matrix
    V = [cos(thetas(i)) sin(thetas(i)); ...
        -sin(thetas(i)) cos(thetas(i))];
    
    % Try to recover the sources using the rotation.
    S_est = V * X_w;
    
    % Calculate the multi-information.
    I(i) =  entropy(S_est(1,:)) + ...
            entropy(S_est(2,:)) - ...
            entropy(S_est);
end

% Plot the multi-information
figure;
hold on; box on
plot(angles, I, '.k-', 'MarkerSize', 16)
xlabel('angle (degrees)');
ylabel('multi information');
xlim([angles(1) angles(end)]);
plot(xlim,zeros(2,1),'--k');

% Plot the original data with the IC's.
figure;
subplot(1,2,1);
hold on; box on;
plot(X(1,:), X(2,:), '.k','MarkerSize', 16);
axis(2*[-2 2 -2 2]); axis square;

% Plot the ICA solution.
[W, S] = ica(X);
subplot(1,2,2);
hold on; box on;
plot(S(1,:), S(2,:), '.k', 'MarkerSize', 16);
axis(2*[-2 2 -2 2]); axis square;

    
    
    
    
   