function [dpri,ccrit] = dprime(pHit,pFA,nTarget,nDistract)

% DPRIME  --  Calculate sensitivity index from signal-detection theory
%
%  USE:
%  dvalue = dprime(pHit,pFA) returns the sensitivity index, using a
%  standard value for correction of rates of 0 and 1.
%  dvalue = dprime(pHit,pFA,nTarget,nDistract) uses the given numbers of
%  targets and distractor trials value for correction.
%  [dvalue,cvalue] = dprime(pHit,pFA) returns also the individual bias
%  criterion.
%
%  Coding by Martin Bckmann-Barthel, Otto-von-Guericke-Universitt Magdeburg
%  The sensitivity index d' measures the distance between the signal and
%  the noise means in standard deviation units. c is the distance of the
%  bias criterion %  from the point where neither response is favored, also
%  in standard units. Positive c values indicate a bias towards high hit
%  rates and false alarm rates
% 
%  Obligatory input variables: 
%  pHit, pFA - hit rate and false alarm rate (max: 1.0, min 0.0)
%  Optional variables: nTarget, nDistract - number of targets and
%  distractor trials (needed for correction of perfect responses)
%  Perfect hit rates (pHit = 1) and FA rates (pFA = 0) are corrected 
%  by -1/(2 nTarget) and +1/(2 nDistract), respectively, if provided
%  cf. Stanislaw H, Todorov N, Behav Res Meth (1999) 31, 137-149, "1/2N rule"
%
%--------------------------------------------------------------------------
%
%  $Revision: 2.0 $  $Date: 2014-07-29 $
%  Handling of perfect response
%  $Revision: 2.01 $  $Date: 2016-08-31 $
%  Help edit
%  $Revision: 3.0 $  $Date: 2017-02-15 $
%  Criterium c (response bias) added
%  $Revision: 3.01 $  $Date: 2017-12-09 $
%  Help edit, input error check

%-- Replace rates equalling zero or one
if nargin < 4 % number of distractor presentations
    nDistract = 1e8; % if not specified, take a very high number
end
if nargin < 3 % number of target presentations
    nTarget = 1e8; % if not specified, take a very high number
end
if pHit > 1 | pFA > 1 
    error('Meaningless probabilities. (Do NOT enter percentage values!)');
end % if
if pHit < 0 | pFA < 0 
    error('Meaningless negative probabilities.');
end % if

if pHit > 0
    pHit = min(pHit,1-.5/nTarget);
else pHit = .5/nTarget;
end % if pHit=0
if pFA < 1
    pFA = max(pFA,.5/nTarget);
else pFA=1-.5/nTarget;
end % if pFA=0

%-- Convert to Z scores, no error checking
zHit = -sqrt(2).*erfcinv(2*pHit);
zFA = -sqrt(2).*erfcinv(2*pFA);
%-- Calculate d-prime
dpri = zHit - zFA ;

%-- Calculate d-prime and bias
dpri = zHit - zFA ;
if nargout > 1
    ccrit = -.5*(zHit + zFA);
end % if nargout
%  Return DPRIME