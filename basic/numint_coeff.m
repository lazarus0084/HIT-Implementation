function XipWip = numint_coeff(N,n,Xd,Xl,a,H)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function organizes the integration points and weights for the
% layered elastic half-space equations

% The integration is always carried out between the first N zeros of the 
% Bessel functions involved. Instead of establishing n integration points,
% only the first interval is integrated using a n-point Gaussian formula.
% The second interval is integrated using a 20-point formula, the third
% interval is integrated using a 10-point formula, and the remaining
% intervals are integrated using a 5-point formula. This formulation is
% chosen to reduce the expense of numerical integration, where most of the
% relevant information is contained between the first two zeros of the 
% Bessel functions. It is recommended to use n >= 30.

% The oranization of integration points and weights are only defined for 
% zeros at B0r locations. This is to further expedite the code performance
% at a later stage (i.e., when inferring arbitrary functions), and have 
% negligeble influence on the resulting response.

% INPUT PARAMETERS
% N  : Number of Bessel roots in integration
% n  : Number of integration points in between each Bessel roots
% Xd : Coordinates of evaluation points (x,y,z)
% Xl : Coordinates of load (x,y)
% a  : Load radii vector
% H  : Depth of layer n-1 out of n layers (i.e. thickness of structure)
%--------------------------------------------------------------------------

% Number of loads (xl) and deformation points (xd)
xl = size(Xl,1); % Number of load points
xd = size(Xd,1); % Number of deformation points

% Alpha - repeated by the number of deformation points all loads influence
alpha = repmat(a/H,[xd,1]);

% Extend number of evluation points
Xl    = repmat(Xl,[xd,1]);  % Extends matrix to number of evaluation points
Xd    = repmat(Xd,[xl,1]);  % Xd is reorganized below to correspond to the
                            % correct Xl content/sequence

dx = repmat(1:xd,xl,1);
dx = dx(:)' + repmat(0:xd:xd*(xl-1),1,xd);
Xd = Xd(dx,:);

% Organize alpha and rho
r = sqrt((Xd(:,1)-Xl(:,1)).^2 + (Xd(:,2)-Xl(:,2)).^2);

% Introduce the parameter rho (notice: since r is a vector, rho is also
% a vector)
rho = r/H;

% Identify roots of Bessel functions
B0r = besselroots(0,N,1); % Bessel 0 roots (B0r)
B1r = besselroots(1,N,1); % Bessel 1 roots (B1r)

% Find nonzero radii/rho (some can be zero if the response at the point of
% loading is considered)
rho_Non0 = find(rho~=0);

n3 = length(rho);      % Number of radii
n4 = length(alpha);    % Number of loads: NB: n3 = n4, due to repmat(..)
                       % commands above

% Organize matrix with integration points:
% We want to know the m0 and m1-values for the roots, i.e. where
% J0(m0*r) = 0, J1(m1*r) = 0 and J1(m1*a) = 0.  These are
% evaluated by taking the above found roots B0r and B1r and as:
% m0 = B0r/r, m1 = B1r/r, m1 = J1/a. The reason why we need the specific m0
% and m1 values is that m-values are included in the the integration as
% separate values. In case r = 0, we have a special case where J0(m0*r) =
% constant and J1(m1*r) = constant. So we do not need parameters for this
% in the integration. We only need integration points for the J1(m1*a).
% So we won't need the m values for the J0(m*r) and J1(m*r) functions. Only
% for J1(m1*a), since it will never be 0. But, below we need the matrix
% size to be consistent. So we put the values to be zero. This does not
% complicate the computations, since (as shown below) we only choose the
% first N nonzero roots in the integration. And we will always have N
% nonzero roots since the radius, a, will always be nonzero. But we still
% need to organize zeros in the roots below, before we can pick N nonzero
% values.
%
% The first xl rows represents the integration points to evaluate the
% deformation of one point due to the xl loads.If a radius is zero,
% division by zero exist

% Define roots matrices with zeros at B0r spaces and B1r roots at the
% remaining spaces
roots_z = [zeros(n3,N) repmat(B1r,[n4,1])./repmat(alpha,[1,N])];

% Insert B0r roots in places where rho > 0 (if rho = 0, B0r/rho = NaN).
if isempty(rho_Non0) ~= 1
    roots_z(rho_Non0,1:N) = repmat(B0r,[length(rho_Non0),1])...
    ./repmat(rho(rho_Non0),[1,N]);
end

% Organize roots in each row in ascending order
roots_z = sort(roots_z')';

% Select the points including up to the first N nonzero zero-value points
% (if we include more than N nonzero zero-value points, we might miss some
% zero points in between <-- because we divide roots by alpha and rho, the
% m-value ranges are different.)

% Define indx vector to select the nonzero roots
indx    = zeros(n3,N); % Organize zero matrix
indx(:) = 1:n3*N;      % Fill up matrix with indeces that assume that no
                       % radius is zero (meaning that roots are present)

% Replace the indeces where the radius is zero with new indeces where
% integration points are nonzero
indx(find(rho==0),:) = indx(find(rho==0),:) + n3*N;

% A special case can appear where alpha = rho, i.e., where the deformation
% point is at the edge of the loading. In this case the
% roots_r = [Br/rho Br/alpha], where Br/rho = Br/alpha. So we have the same
% integration points appearing twice. These are arranged side by side after
% we have used the 'sort'-function. So in this case we have to replace the
% index with another index to ensure we do not take out the same values
% twice. This is done below where we add a vector to each of these rows
% given as [0 1 2 3 .... N-2 N-1]*n3, i.e., the collums are shifted (we
% pick out every second column instead of the first N columns).

indx(find(rho==alpha),:) = indx(find(rho==alpha),:) + (0:N-1)*n3;

% % if bug in above lone, use this
% if length(find(rho==alpha)>0)
%     indx(find(rho==alpha),:) = indx(find(rho==alpha),:) + (0:N-1)*n3;
% end

% Pick out the nonzero roots
seq_z = [zeros(n3,1) roots_z(indx)];

% Distribution of n integration points and weights of the first interval 
[xip_z1 , wip_z1]=lookup_gauss(n,seq_z(:,1),seq_z(:,2));   

% Distribution of 20 integration points and weights of the second interval
[xip_z2 , wip_z2]=lookup_gauss(20,seq_z(:,2),seq_z(:,3));   

% Distribution of 10 integration points and weights of the third interval
[xip_z3 , wip_z3]=lookup_gauss(10,seq_z(:,3),seq_z(:,4));   

% Define remaining besselzeroes intervals
seq_z1 = seq_z(:,4:end-1)'; % <-- Note: transpose, in order for each collumn
                            % to represent the integration points for ONE
                            % deformation due to ONE load
seq_z2 = seq_z(:,5:end)';

% Distribution of 5 integration points and weights of the remaining intervals
[xip_z , wip_z]=lookup_gauss(5,seq_z1(:),seq_z2(:)); % z: ip for z-direction

% Reorganize integration points and weights - Every row
% (after reorganization) now refers to a single radius/point-to-load case
xip_z = xip_z';                        % First we transform the matrices
xip_z = reshape(xip_z(:),[],n3)'; 
wip_z = wip_z';
wip_z = reshape(wip_z(:),[],n3)';

% Collect all integration points
xip_z = [xip_z1 xip_z2 xip_z3 xip_z];
xip_r = xip_z;                         % Uses only the integration point 
                                       % defined for zeros at B0r locations

% Collect all weights
wip_z = [wip_z1 wip_z2 wip_z3 wip_z];
wip_r = wip_z;

XipWip = [xip_z; wip_z; xip_r; wip_r];
