% Demo simulation of DGCMG dynamics

clear
clc

q2  = 0.3;
q3  = -0.2;
q2d = 0.5;
q3d = 0.4;
q4d = 100;

[M,C,W] = computeGyroDynamics(q2,q3,q2d,q3d,q4d);

disp("Inertia Matrix M:")
disp(M)

disp("Coriolis Matrix C:")
disp(C)

disp("Inverse Inertia Matrix W:")
disp(W)
