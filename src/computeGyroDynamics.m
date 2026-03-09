function [M, C, W] = computeGyroDynamics(q2, q3, q2d, q3d, q4d)
%COMPUTEGYRODYNAMICS Compute inertia, Coriolis, and inverse inertia matrices
% for a 4-DOF double-gimbal control moment gyroscope (DGCMG) model.
%
% Inputs
%   q2   - Gimbal angle 2 [rad]
%   q3   - Gimbal angle 3 [rad]
%   q2d  - Gimbal angular velocity 2 [rad/s]
%   q3d  - Gimbal angular velocity 3 [rad/s]
%   q4d  - Rotor angular velocity 4 [rad/s]
%
% Outputs
%   M    - 4x4 generalized inertia matrix
%   C    - 4x4 Coriolis/centrifugal matrix
%   W    - 4x4 inverse inertia matrix, inv(M)
%
% Notes
%   - Beta is treated as a constant rotor spin rate.
%   - The Coriolis matrix is assembled from Christoffel-like coefficient
%     matrices L1-L4.
%   - This implementation is written for readability, traceability,
%     and integration into simulation / control pipelines.

    %-----------------------------%
    % Input validation
    %-----------------------------%
    validateattributes(q2,  {'numeric'}, {'real','scalar','finite'}, mfilename, 'q2',  1);
    validateattributes(q3,  {'numeric'}, {'real','scalar','finite'}, mfilename, 'q3',  2);
    validateattributes(q2d, {'numeric'}, {'real','scalar','finite'}, mfilename, 'q2d', 3);
    validateattributes(q3d, {'numeric'}, {'real','scalar','finite'}, mfilename, 'q3d', 4);
    validateattributes(q4d, {'numeric'}, {'real','scalar','finite'}, mfilename, 'q4d', 5);

    %-----------------------------%
    % Physical parameters
    %-----------------------------%
    p = getModelParameters();

    % Constant spin speed [rad/s]
    beta = 750 * pi / 30;

    %-----------------------------%
    % Precompute reusable terms
    %-----------------------------%
    s2 = sin(q2);
    c2 = cos(q2);
    s3 = sin(q3);
    c3 = cos(q3);

    s2_sq = s2^2;
    c2_sq = c2^2;
    s3_sq = s3^2;
    c3_sq = c3^2;

    s2c2 = s2 * c2;

    % Lumped coefficients
    a1 = p.Jc - p.Kc;
    a2 = p.Jd - p.Id;
    a3 = p.Id - p.Jc - p.Jd + p.Kc;
    a4 = p.Ic + p.Id;
    a5 = p.Ib + p.Ic - p.Kb - p.Kc;

    %-----------------------------%
    % Inertia matrix M(q)
    %-----------------------------%
    Ma = [0, 0, 0, 0;
          0, 0, 0, 0;
          0, 0, 0, 0;
          0, 0, 0, p.Ka];

    Mb = [0, 0, 0, 0;
          0, 0, 0, 0;
          0, 0, p.Jb, 0;
          0, 0, 0, p.Ib*s3_sq + p.Kb*c3_sq];

    Mc = [0, 0, 0, 0;
          0, p.Ic, 0, -p.Ic*s3;
          0, 0, p.Jc*c2_sq + p.Kc*s2_sq, a1*s2c2*c3;
          0, -p.Ic*s3, a1*s2c2*c3, p.Ic*s3_sq + (p.Jc*s2_sq + p.Kc*c2_sq)*c3_sq];

    Md = [p.Jd, 0, p.Jd*c2, p.Jd*s2*c3;
          0, p.Id, 0, -p.Id*s3;
          p.Jd*c2, 0, p.Id*s2_sq + p.Jd*c2_sq, a2*s2c2*c3;
          p.Jd*s2*c3, -p.Id*s3, a2*s2c2*c3, p.Id*s3_sq + c3_sq*(p.Id*c2_sq + p.Jd*s2_sq)];

    M = Ma + Mb + Mc + Md;

    %-----------------------------%
    % Christoffel-like coefficient matrices
    %-----------------------------%
    L1 = 0.5 * [0, 0, 0, 0;
                0, 0, -p.Jd*s2, p.Jd*c2*c3;
                0, -p.Jd*s2, 0, -p.Jd*s2*s3;
                0, p.Jd*c2*c3, -p.Jd*s2*s3, 0];

    L2 = 0.5 * [0, 0, p.Jd*s2, -p.Jd*c2*c3;
                0, 0, 0, 0;
                p.Jd*s2, 0, -2*a3*s2c2, a3*(c2_sq*c3 - s2_sq*c3) - a4*c3;
                -p.Jd*c2*c3, 0, a3*(c2_sq*c3 - s2_sq*c3) - a4*c3, 2*a3*c2*c3_sq*s2];

    L3 = 0.5 * [0, -p.Jd*s2, 0, p.Jd*s2*s3;
                -p.Jd*s2, 0, 2*a3*s2c2, a4*c3 + a3*(c3*s2_sq - c2_sq*c3);
                0, 2*a3*s2c2, 0, 0;
                p.Jd*s2*s3, a4*c3 + a3*(c3*s2_sq - c2_sq*c3), 0, -2*(a5 + a3*s2_sq)*c3*s3];

    L4 = 0.5 * [0, p.Jd*c2*c3, -p.Jd*s2*s3, 0;
                p.Jd*c2*c3, 0, a3*(c3*s2_sq - c2_sq*c3) - a4*c3, -2*a3*c2*c3_sq*s2;
                -p.Jd*s2*s3, a3*(c3*s2_sq - c2_sq*c3) - a4*c3, 2*a3*c2*s2*s3, 2*(a5 + a3*s2_sq)*c3*s3;
                0, -2*a3*c2*c3_sq*s2, 2*(a5 + a3*s2_sq)*c3*s3, 0];

    %-----------------------------%
    % Generalized velocity vector
    %-----------------------------%
    qdRow = [beta, q2d, q3d, q4d];
    qdCol = qdRow.';

    %-----------------------------%
    % Coriolis / centrifugal matrix C(q, qd)
    %-----------------------------%
    velocityBlock = blkdiag(qdRow, qdRow, qdRow, qdRow);
    Lstack = [L1; L2; L3; L4];

    C = velocityBlock * Lstack;

    %-----------------------------%
    % Inverse inertia matrix W = inv(M)
    % Prefer matrix division for numerical robustness
    %-----------------------------%
    W = M \ eye(4);

    % Optional numerical health checks
    if any(~isfinite(M), 'all') || any(~isfinite(C), 'all') || any(~isfinite(W), 'all')
        error('computeGyroDynamics:NonFiniteResult', ...
              'Non-finite value detected in computed system matrices.');
    end

    % Optional consistency check
    if rcond(M) < 1e-12
        warning('computeGyroDynamics:IllConditionedMassMatrix', ...
                'Mass matrix M is close to singular or poorly conditioned.');
    end
end


function p = getModelParameters()
%GETMODELPARAMETERS Return fixed physical parameters for the DGCMG model.

    p.Ka = 0.0375;
    p.Kb = 0.0213;
    p.Kc = 0.0027;

    p.Ib = 0.0037;
    p.Ic = 0.0010;
    p.Id = 0.0028;

    p.Jb = 0.0179;
    p.Jc = 0.0018;
    p.Jd = 0.0056;
end
