%<removed> function [xdot,U,M] = remus100(x,ui,Vc,betaVc,w_c)
function [xdot,U,M,tau_out] = remus100_consistent(x,ui,Vc,betaVc,w_c,opt)  %<added>
% remus100_consistent.m: modified copy of MSS CRAFT/AUV/models/remus100.m (T. I. Fossen),  %<added>
% MSS revision cc07579 (https://github.com/cybergalactic/MSS).  %<added>
% Change: one value per physical quantity (one water density, one geometry,  %<added>
% the mass given), each changed line commented 'one value:':  %<added>
%   - one water density rho (line 98, 1026 kg/m^3) reaches every call: the  %<added>
%     added mass (imlay61.m 31 holds 1026), the hull lift and drag  %<added>
%     (forceLiftDrag.m 26, 1026) and the cross-flow (crossFlowDrag.m 36, 1025)  %<added>
%     through marked copies that take rho as an input;  %<added>
%   - one geometry, the hull's L_auv = 1.6 m and D_auv = 0.19 m (lines 131-132),  %<added>
%     for the mass, the added mass, CD_0 and the fin positions too (lines  %<added>
%     134-135 scale both by 1.0096 for the mass only);  %<added>
%   - the mass given, 31.9 kg (line 3), not 1025 kg/m^3 times the scaled  %<added>
%     spheroid volume (spheroid.m 35-36); W = m g_mu and B = W as line 214.  %<added>
% Options of this copy (6th argument, a struct; defaults give the one-value model):  %<added>
%   opt.rho       the one water density (default line 98)  %<added>
%   opt.mass      the given mass (default 31.9 kg, line 3)  %<added>
%   opt.munk      true keeps the added-mass Coriolis couplings of lines 207-210  %<added>
%   opt.tau_ext   6x1 external generalized force (N, N m; BODY, CO) in place of  %<added>
%                 the actuator wrench: the bare hull (ui must then be zeros(3,1))  %<added>
%   opt.geometry  'one' (default) or 'mss': MSS's scaled semi-axes and the mass  %<added>
%                 of the displaced water at rho (spheroid.m 36), i.e. MSS with  %<added>
%                 one density only (history: separates density from geometry)  %<added>
% Fourth output tau_out: the actuator wrench of lines 248-254.  %<added>
% Every changed line is marked: a line starting with '%<removed> ' is the  %<added>
% original line commented out; a line ending with '%<added>' is new. Deleting  %<added>
% the added lines and un-commenting the removed ones gives the original file  %<added>
% byte for byte (checked by generate_remus100_consistent.m before it runs).  %<added>
% Modified by: Enio Krizman, 2026-10-07. The original file's licence:  %<added>
%  %<added>
% MIT License  %<added>
%  %<added>
% Copyright (c) 2004 Thor I. Fossen  %<added>
%  %<added>
% Permission is hereby granted, free of charge, to any person obtaining a copy  %<added>
% of this software and associated documentation files (the "Software"), to deal  %<added>
% in the Software without restriction, including without limitation the rights  %<added>
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell  %<added>
% copies of the Software, and to permit persons to whom the Software is  %<added>
% furnished to do so, subject to the following conditions:  %<added>
%  %<added>
% The above copyright notice and this permission notice shall be included in all  %<added>
% copies or substantial portions of the Software.  %<added>
%  %<added>
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  %<added>
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,  %<added>
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE  %<added>
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER  %<added>
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  %<added>
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  %<added>
% SOFTWARE.  %<added>
%  %<added>
% The length of the Remus 100 AUV is 1.6 m, the cylinder diameter is 19 cm  
% and the mass of the vehicle is 31.9 kg. The maximum speed of 2.5 m/s is 
% obtained when the propeller runs at 1525 rpm in zero currents. The
% function returns the time derivative xdot of the state vector: 
%
%   x = [ u v w p q r x y z phi theta psi ]', alternatively 
%   x = [ u v w p q r x y z eta eps1 eps2 eps3 ]' 
%
% in addition to the speed U in m/s (optionally). The state vector can be 
% of dimension 12 (Euler angles) or 13 (unit quaternions):
%
%   u:       Surge velocity          (m/s)
%   v:       Sway velocity           (m/s)
%   w:       Heave velocity          (m/s)
%   p:       Roll rate               (rad/s)
%   q:       Pitch rate              (rad/s)
%   r:       Yaw rate                (rad/s)
%   x:       North position          (m)
%   y:       East position           (m)
%   z:       Downwards position      (m)
%   phi:     Roll angle              (rad)       
%   theta:   Pitch angle             (rad)
%   psi:     Yaw angle               (rad)
% 
% For the unit quaternion representation, the last three arguments of the 
% x-vector, the Euler angles (phi, theta, psi), are replaced by the unit 
% quaternion quat = [eta, eps1, eps2, eps3]'. This increases the dimension of
% the state vector from 12 to 13.
%
% The control inputs are one tail rudder, two stern planes and a single-screw 
% propeller:
%
%   ui = [ delta_r delta_s n ]'  where
%
%    delta_r:   Rudder angle (rad)
%    delta_s:   Stern plane angle (rad) 
%    n_p:       Propeller revolution (RPM)
%
% The arguments Vc (m/s), betaVc (rad), w_c (m/s) are optional arguments for 
% ocean currents
%
%    v_c = [ Vc * cos(betaVc - psi), Vc * sin( betaVc - psi), w_c ]  
% 
% Example usage: 
%   [~,~,M] = remus100()                            : Return the 6x6 mass matrix M
%   [xdot,U] = remus100(x,ui,Vc,betaVc,alphaVc,w_c) : 3-D ocean currents
%   [xdot,U] = remus100(x,ui,Vc,betaVc,alphaVc)     : 2-D ocean currents
%   [xdot,U] = remus100(x,ui)                       : No ocean currents
%   xdot = remus100(x,ui)                           : No ocean currents
%
% Author:    Thor I. Fossen
% Date:      2021-05-27
% Revisions:
%   2021-08-24  Ocean currents are now expressed in NED
%   2021-10-21  imlay61.m is called using the relative velocity
%   2021-12-30  Added the time derivative of the current velocity
%   2022-02-01  Updated lift and drag forces
%   2022-05-06  Calibration of drag and propulsion forces using data from 
%               Allen et al. (2000)
%   2022-06-08  Added compatibility for unit quaternions in addition to the 
%               Euler angle representation
%   2022-10-16  Added vertical currents
%   2023-05-02  Corrected the rudder area A_r
%   2023-10-07  Scaled down the propeller roll-induced moment
%   2024-02-09  Updated rudder and stern-plane areas
%   2024-02-13  Calibration of the model parameters
%   2024-06-23  Corrected sign of tau(5). A positive delta_s will result in a
%               negative pitch and negative Z_s. Hence, tau(5) = -x_s * Z_s 
%               when x_s < 0.
%   2025-04-25 Added empty call: [~,~,M] = remus100(), and minor bug fixes.
%   2025-06-09 Change cross-flow drag to cylinder model (M. Seidl).
%   2026-08-26 Use through-water speed for propulsion and damping (J. Harvey).
%
% References: 
%   B. Allen, W. S. Vorus and T. Prestero, "Propulsion system 
%       performance enhancements on REMUS AUVs," OCEANS 2000 MTS/IEEE 
%       Conference and Exhibition. Conference Proceedings, 2000, 
%       pp. 1869-1873 vol.3, doi: 10.1109/OCEANS.2000.882209.
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and
%       Motion Control. 3rd. Edition, Wiley. URL: www.fossen.biz/wiley   

if nargin == 0  % [~,~,M] = remus100()
    x = zeros(12,1); ui = zeros(3,1); Vc = 0; betaVc = 0; w_c = 0;
end

if (nargin == 2), Vc = 0; betaVc = 0; w_c = 0; end % No ocean currents
if (nargin == 4), w_c = 0; end % No vertical ocean currents
if nargin < 6, opt = struct(); end  % one value: options of this copy  %<added>
if ~isfield(opt, 'geometry'), opt.geometry = 'one'; end  % one value: one geometry by default  %<added>
if ~isfield(opt, 'munk'), opt.munk = false; end  % one value: MSS's zeroed couplings by default  %<added>

if (length(ui) ~= 3),error('u-vector must have dimension 3!'); end
if (length(x) ~= 12 && length(x) ~= 13)
    error('x-vector must have dimension 12 or 13'); 
end

% Constants
mu = deg2rad(63.446827); % Latitude for Trondheim, Norway (deg)
g_mu = gravity(mu);      % Gravity vector (m/s2)
rho = 1026;              % Density of water (m/s2)
if isfield(opt, 'rho'), rho = opt.rho; end  % one value: the one water density, passed to every call below  %<added>

% 6x1 velocity vector, yaw angle and control inputs
nu = x(1:6); 
if length(x) == 12 % Euler angles
    psi = x(12);
else % Convert unit quaternion to yaw angle
    psi = atan2(2*(x(10)*x(13) + x(11)*x(12)), 1 - 2*(x(12)^2 + x(13)^2));
end

% Amplitude saturation of the control signals
delta_max = deg2rad(20); % Maximum rudder and stern angles (rad)
n_max = 1525;            % Maximum propeller speed (RPM)

% Amplitude saturation of the control signals
delta_r = sat(ui(1), delta_max);    % Saturated tail rudder (rad)
delta_s = sat(ui(2), delta_max);    % Saturated Stern plane (rad)
n_p = sat(ui(3),n_max) / 60;        % Saturated propeller speed (rps)

% Ocean currents expressed in BODY
u_c = Vc * cos( betaVc - psi );                               
v_c = Vc * sin( betaVc - psi );   

nu_c = [u_c v_c w_c 0 0 0]'; % Ocean current velocities
Dnu_c = [nu(6)*v_c -nu(6)*u_c 0 0 0 0]'; % Time derivative of nu_c

% Relative velocities/speed, angle of attack and vehicle speed
nu_r = nu - nu_c;                                 % Relative velocity
alpha = atan2( nu_r(3), nu_r(1) );                % Angle of attack (rad)
U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );  % Relative speed (m/s)
U  = sqrt( nu(1)^2 + nu(2)^2 + nu(3)^2 );         % Speed (m/s)

% AUV model parameters; Fossen (2027, Chapter 8) and Allen et al. (2000)
L_auv = 1.6;             % AUV length (m)
D_auv = 0.19;            % AUV diamater (m)
S = 0.7 * L_auv * D_auv; % Planform area S = 70% of rectangle L_auv * D_auv
a = 1.0096 * L_auv/2;    % Scaled spheroid semi-axes a and b to obtain m = 31.9 kg
b = 1.0096 * D_auv/2;                   
if strcmp(opt.geometry, 'one'), a = L_auv/2; b = D_auv/2; end  % one value: the hull's own L_auv, D_auv (lines 131-132)  %<added>
m_given = 31.9;  % one value: the mass, given (line 3: 'the mass of the vehicle is 31.9 kg')  %<added>
if isfield(opt, 'mass'), m_given = opt.mass; end  % one value: the given mass as an option  %<added>
if strcmp(opt.geometry, 'mss'), m_given = 4/3 * pi * rho * a * b^2; end  % one value: MSS geometry, displaced mass at the one rho (spheroid.m 36)  %<added>
r44 = 0.3;               % Added moment of inertia in roll: A44 = r44 * Ix
r_bG = [ 0 0 0.02 ]';    % CG w.r.t. to the CO
r_bB = [ 0 0 0 ]';       % CB w.r.t. to the CO

% Parasitic drag coefficient CD_0, i.e. zero lift and alpha = 0
% F_drag = 0.5 * rho * Cd * (pi * b^2)   
% F_drag = 0.5 * rho * CD_0 * S
Cd = 0.42;                              % From Allen et al. (2000)
CD_0 = Cd * pi * b^2 / S;

% Propeller coeffs. KT and KQ are computed as a function of advance no.
% Ja = Va/(n*D_prop) where Va = (1-w)*U_r = 0.944 * U_r; Allen et al. (2000)
D_prop = 0.14;   % Propeller diameter corresponding to 5.5 inches
t_prop = 0.1;    % Thrust deduction number
Va = 0.944 * U_r;  % Advance speed (m/s)

% Ja_max = 0.944 * 2.5 / (0.14 * 1525/60) = 0.6632
Ja_max = 0.6632;
        
% Single-screw propeller with 3 blades and blade-area ratio = 0.718.    
% >> [KT_0, KQ_0] = wageningen(0,1,0.718,3)
KT_0 = 0.4566;
KQ_0 = 0.0700;
% >> [KT_max, KQ_max] = wageningen(0.6632,1,0.718,3) 
KT_max = 0.1798;
KQ_max = 0.0312;
        
% Propeller thrust and propeller-induced roll moment
% Linear approximations for positive Ja values
% KT ~= KT_0 + (KT_max-KT_0)/Ja_max * Ja   
% KQ ~= KQ_0 + (KQ_max-KQ_0)/Ja_max * Ja        
if n_p > 0   
    % Forward thrust
    X_prop = rho * D_prop^4 * ( ... 
        KT_0 * abs(n_p) * n_p + (KT_max-KT_0)/Ja_max * (Va/D_prop) * abs(n_p) );        
    K_prop = rho * D_prop^5 * ( ...
        KQ_0 * abs(n_p) * n_p + (KQ_max-KQ_0)/Ja_max * (Va/D_prop) * abs(n_p) );                 
else      
    % Reverse thrust (braking)   
    X_prop = rho * D_prop^4 * KT_0 * abs(n_p) * n_p; 
    K_prop = rho * D_prop^5 * KQ_0 * abs(n_p) * n_p;            
end            

S_fin = 0.00665;         % Fin area

% Tail rudder
CL_delta_r = 0.5;        % Rudder lift coefficient (-)
A_r = 2 * S_fin;         % Rudder area (m2)
x_r = -a;                % Rudder x-position (m)

% Stern plane (double)
CL_delta_s = 0.7;        % Stern-plane lift coefficient (-)
A_s = 2 * S_fin;         % Stern-plane area (m2)
x_s = -a;                % Stern-plane z-position (m)

% Low-speed linear damping matrix parameters
T1 = 20;                 % Time constant in surge (s)
T2 = 20;                 % Time constant in sway (s)
zeta4 = 0.3;             % Relative damping ratio in roll
zeta5 = 0.8;             % Relative damping ratio in pitch
T6 = 1;                  % Time constant in yaw (s)

% Rigid-body mass and hydrodynamic added mass
%<removed> [MRB,CRB] = spheroid(a,b,nu(4:6),r_bG);
[MRB,CRB] = spheroid_consistent(a,b,nu(4:6),r_bG,m_given);  % one value: the given mass, not rho = 1025 of spheroid.m 35  %<added>
%<removed> [MA,CA] = imlay61(a, b, nu_r, r44);
[MA,CA] = imlay61_consistent(a, b, nu_r, r44, rho);  % one value: rho, not 1026 of imlay61.m 31  %<added>

% The reduced-order REMUS model does not retain a complete set of measured
% maneuvering derivatives for the hull, fins, and appendages. The selected
% added-mass Coriolis couplings are therefore omitted to avoid retaining an
% unbalanced Munk-moment model. This is a modeling assumption; the
% corresponding physical terms are not generally zero for a bare body.
if ~opt.munk  % one value: opt.munk keeps C_A = m2c(M_A, nu_r) in full  %<added>
CA(5,3) = 0; CA(3,5) = 0;  % Heave-pitch coupling
CA(5,1) = 0; CA(1,5) = 0;  % Surge-pitch Munk coupling
CA(6,1) = 0; CA(1,6) = 0;  % Yaw-related Munk couplings
CA(6,2) = 0; CA(2,6) = 0;
end  % one value: end of the option  %<added>

M = MRB + MA;
C = CRB + CA;
m = MRB(1,1); W = m * g_mu; B = W;

% Dissipative forces and moments
D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bG' r_bB']);
D(1,1) = D(1,1) * exp(-3 * U_r); % Vanish at high through-water speed

%<removed> tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);
tau_liftdrag = forceLiftDrag_consistent(D_auv,S,CD_0,alpha,U_r,rho);  % one value: rho, not 1026 of forceLiftDrag.m 26  %<added>
%<removed> tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r,'cylinder');
tau_crossflow = crossFlowDrag_consistent(L_auv,D_auv,D_auv,nu_r,'cylinder',rho);  % one value: rho, not 1025 of crossFlowDrag.m 36  %<added>

% Kinematics
if (length(x) == 13)
    [J,R] = quatern(x(10:13));
else
    [J,R] = eulerang(x(10),x(11),x(12));
end

% Restoring forces and moments
g = gRvect(W,B,R,r_bG,r_bB);

% Horizontal- and vertical-plane relative speed
U_rh = sqrt( nu_r(1)^2 + nu_r(2)^2 );  
U_rv = sqrt( nu_r(1)^2 + nu_r(3)^2 );  

% Rudder and stern-plane drag
X_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r^2; 
X_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s^2;

% Rudder sway force 
Y_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r;

% Stern-plane heave force
Z_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s;

% Generalized propulsion force vector
tau = zeros(6,1);                                
tau(1) = (1-t_prop) * X_prop + X_r + X_s;
tau(2) = Y_r;
tau(3) = Z_s;
tau(4) = K_prop / 10; % Scaled down by a factor of 10 to match exp. results
tau(5) = -x_s * Z_s;  
tau(6) = x_r * Y_r;
tau_out = tau;  % one value: the actuator wrench, returned  %<added>
if isfield(opt, 'tau_ext'), tau = opt.tau_ext(:); end  % one value: bare hull, an external wrench in place of the actuators  %<added>

% State-space model
xdot = [ Dnu_c + M \ ...
            (tau + tau_liftdrag + tau_crossflow - C * nu_r - D * nu_r  - g)
         J * nu ];
