%<removed> function tau_crossflow = crossFlowDrag(L,B,T,nu_r,drag_model)
function tau_crossflow = crossFlowDrag_consistent(L,B,T,nu_r,drag_model,rho)  %<added>
% crossFlowDrag_consistent.m: modified copy of MSS LIBRARY/modeling/crossFlowDrag.m (T. I. Fossen),  %<added>
% MSS revision cc07579 (https://github.com/cybergalactic/MSS).  %<added>
% Change: the water density rho is an input, not 1025 kg/m^3 (line 36).  %<added>
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
% tau_crossflow = crossFlowDrag(L,B,T,nu_r,drag_model) computes the cross-flow 
% drag integrals for a marine craft using strip theory. Application:
%
%    M * d/dt nu_r + C(nu_r) * nu_r + D * nu_r + g(eta) = tau + tau_crossflow
%
% Inputs: 
%   L:  length (m)
%   B:  beam (m)
%   T:  draft (m)
%   nu_r = [u-u_c, v-v_c, w-w_c, p, q, r]': Relative velocity vector
%   drag_model (Optionally) : 'Hoerner' (default), 'cylinder'
%
% Output: 
%   tau_crossflow = [0 Yh Zh 0 Mh Nh]: 6-DOF cross-flow drag vector
%
% References:
%   Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and Motion
%   Control, 3rd ed., John Wiley & Sons Ltd., Chichester, UK.
%
%   DNV (2025). DNV-RP-C205: Environmental Conditions and Environmental
%   Loads. Recommended Practice, Edition 2025-04, amended 2026-03.
%
% Author:     Thor I. Fossen 
% Date:       25 Apr 2021
% Revisions:  
%   2021-01-30 : Extended to include heave and pitch for AUVs
%   2025-06-09 : Make different drag models selectable (M. Seidl)
%   2026-08-26 : Evaluate 20 strips at their midpoints (J. Harvey)

% Hoerner is the default drag model for monohull
if nargin == 4
    drag_model = 'Hoerner'; 
end 

%<removed> rho = 1025;         % Density of water (kg/m3)
nStrips = 20;       % Number of hull strips
dx = L / nStrips;   % Strip width (m)

switch drag_model

    case 'Hoerner'
        % 2-D drag coefficient based on Hoerner's curve for monohull
        Cd_2D = Hoerner(B, T); 

    case 'cylinder'
        % 2D drag coefficient based on cylinder data from DNV-RP-C205
        Cd_2D = cylinderDrag(L, B, nu_r); 

    otherwise
        error('Unsupported drag model %s.', drag_model)
end

Yh = 0; Zh = 0; Mh = 0; Nh = 0;
for i = 1:nStrips
    xL = -L/2 + (i - 0.5) * dx;                       % Strip midpoint
    v_r = nu_r(2);                                    % Relative sway velocity
    w_r = nu_r(3);                                    % Relative heave velocity
    q = nu_r(5);                                      % Pitch rate
    r = nu_r(6);                                      % Yaw rate
    U_h = abs(v_r + xL * r) * (v_r + xL * r);
    U_v = abs(w_r + xL * q) * (w_r + xL * q);   
    Yh = Yh - 0.5 * rho * T * Cd_2D * U_h * dx;       % Sway force
    Zh = Zh - 0.5 * rho * T * Cd_2D * U_v * dx;       % Heave force  
    Mh = Mh - 0.5 * rho * T * Cd_2D * xL * U_v * dx;  % Pitch moment    
    Nh = Nh - 0.5 * rho * T * Cd_2D * xL * U_h * dx;  % Yaw moment
end

tau_crossflow = [0 Yh Zh 0 Mh Nh]';

end
