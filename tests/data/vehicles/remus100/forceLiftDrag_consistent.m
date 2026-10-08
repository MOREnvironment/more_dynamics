%<removed> function tau_liftdrag = forceLiftDrag(b,S,CD_0,alpha,U_r)
function tau_liftdrag = forceLiftDrag_consistent(b,S,CD_0,alpha,U_r,rho)  %<added>
% forceLiftDrag_consistent.m: modified copy of MSS LIBRARY/modeling/forceLiftDrag.m (T. I. Fossen),  %<added>
% MSS revision cc07579 (https://github.com/cybergalactic/MSS).  %<added>
% Change: the water density rho is an input, not 1026 kg/m^3 (line 26).  %<added>
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
% tau_liftdrag = forceLiftDrag(b,S,CD_0,alpha,Ur) computes the hydrodynamic
% lift and drag forces of a submerged "wing profile" for varying angle of
% attack (Beard and McLain 2012). Application:
%
%  M d/dt nu_r + C(nu_r)*nu_r + D*nu_r + g(eta) = tau + tau_liftdrag
%
% Output:
%  tau_liftdrag:  6x1 generalized force vector
%
% Inputs:
%  b:       wing span (m)
%  S:       wing area (m^2)
%  CD_0:    parasitic drag (alpha = 0), typically 0.1-0.2 for a streamlined body
%  alpha:   angle of attack, scalar or vector (rad)
%  U_r:     relative speed (m/s)
%
% Example:
%
% Cylinder-shaped AUV with length L = 1.8, diameter D = 0.2 and CD_0 = 0.1:
%    tau_liftdrag = forceLiftDrag(0.2, 1.8*0.2, 0.1, alpha, U_r)
% 
% Author:    Thor I. Fossen
% Date:      25 April 2021 

%<removed> rho = 1026;

[CL,CD] = coeffLiftDrag(b,S,CD_0,alpha,0);

F_drag = 1/2 * rho * U_r^2 * S * CD;    % drag force
F_lift = 1/2 * rho * U_r^2 * S * CL;    % lift force

% transform from FLOW axes to BODY axes using angle of attack
tau_liftdrag = [...
    cos(alpha) * (-F_drag) - sin(alpha) * (-F_lift)
    0
    sin(alpha) * (-F_drag) + cos(alpha) * (-F_lift)
    0
    0
    0 ];