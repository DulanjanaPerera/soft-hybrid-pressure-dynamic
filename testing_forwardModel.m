% Testing the forward and inverse pressure models
% 
% 
% 

clear

A = pi*(0.013/2)^2;
r = 0.013;
K = 0.0155;
L = 0.2778;

p=[0.0001;0.01;0.01];
phi = linspace(0.0001, pi/3, 500);
theta = linspace(0, 0, 500);
bend = zeros(size(phi));
dir = zeros(size(phi));
p_mat = zeros(3,length(phi));
p_mat(:,1) = p;

frw = zeros(2, length(phi));

temp_c = pressure2config(p, r, K, A);
theta(1) = temp_c(1);
phi(1) = wrapToPi(temp_c(2));


frw(:,1) = [theta(1); phi(1)];

for i=1:(length(phi))
    
    temp_c = pressure2config(p, r, K, A);
    dir(i) = temp_c(1);
    bend(i) = wrapToPi(temp_c(2)); 

    if (i ~= length(phi))
        J = pressure2configJacobian(p, A, r, K);
        J_inv = pinv(J);
        delta_p = J \ ([theta(i+1); phi(i+1)]-[dir(i); bend(i)]);
        p = p_mat(:, i) + [-(delta_p(1) + delta_p(2)); delta_p(1); delta_p(2)];
        % p = PositivePressures(p);
        p_mat(:, i+1) = p; 
    
        % forward model
        frw(:,i+1) = wrapToPi(frw(:,i) + J*(p_mat(2:3, i+1) - p_mat(2:3, i)));
    end
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
plot((rad2deg(dir)), '*')
hold on
% plot(rad2deg(frw(1,:)), '+')
plot((rad2deg(theta)),'--')
hold off
grid on
axis tight
legend theta forward 'theta desired'

figure(3)
plot((rad2deg(bend)), '*')
hold on
% plot(rad2deg(frw(2,:)), '+')
plot((rad2deg(phi)), '--')
hold off
grid on
axis tight
legend phi forward 'phi desired'

figure(4)
plot(wrapTo360(rad2deg(wrapTo2Pi(abs(dir - theta)))))
hold off
grid on
axis tight
legend 'theta difference'