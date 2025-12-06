% Testing the forward and inverse pressure models

clear

A = pi*(0.013/2)^2;
r = 0.013;
K = 0.0226;
L = 0.2778;

p=[0.0001;0.01;0.01];
phi = linspace(0.1, 0.1, 500);
t = linspace(0,1,500);
theta = wrapToPi(pi-pi*sawtooth(2*pi*t));
% theta = linspace(0, 0, 500);
p_mat = zeros(3,length(phi));
p_mat(:,1) = p;

frw = zeros(2, length(phi));

temp_c = pressure2config(p, r, K, A);
theta(1) = temp_c(1);
phi(1) = wrapToPi(temp_c(2));


frw(:,1) = [theta(1); phi(1)];

for i=1:(length(phi))
    
    [p_est, error_c] = config2pressure(theta(i), phi(i), r, K, A);
    p_mat(:, i) = p_est;

    % forward model
    frw(:,i) = pressure2config(p_est, r, K, A);

end

figure(1);
plot((1:1:length(phi)),p_mat(1,:),'*');
hold on
plot((1:1:length(phi)),p_mat(2,:),'o', 'Color','k');
plot((1:1:length(phi)),p_mat(3,:),'--');
hold off
grid on
axis tight
legend p1 p2 p3

figure(2)
plot((rad2deg(frw(1,:))), '*')
hold on
% plot(rad2deg(frw(1,:)), '+')
plot((rad2deg(theta)),'--')
hold off
grid on
axis tight
legend 'theta forward' 'theta desired'

figure(3)
plot((rad2deg(frw(2,:))), '*')
hold on
% plot(rad2deg(frw(2,:)), '+')
plot((rad2deg(phi)), '--')
hold off
grid on
axis tight
legend 'phi forward' 'phi desired'

figure(4)
plot(wrapTo360(rad2deg(wrapTo2Pi(abs(frw(1,:) - theta)))))
hold off
grid on
axis tight
legend 'theta difference'