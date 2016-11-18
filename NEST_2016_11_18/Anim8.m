load M_vs_t

% Data stored:
% save M_vs_t M_x M_y M_z M_x_tot M_x_tot M_x_tot Happ_x Happ_y Happ_z t_solved x y z

Ms=max(sqrt(M_x(1,:).*M_x(1,:)+M_y(1,:).*M_y(1,:)+M_z(1,:).*M_z(1,:)))
delta=min(abs(diff(x)));
figure(10)
clf
M_x=M_x/Ms;
M_y=M_y/Ms;
M_z=M_z/Ms;
for i=2500:50:5000
    clf
    quiver3 (x,y,z,M_x(i,:)',M_y(i,:)',M_z(i,:)')
    axis equal
    pause
end
