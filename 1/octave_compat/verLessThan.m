function tf = verLessThan(pkg, ver) %#ok<INUSD>
%VERLESSTHAN  GNU Octave shim for redbKIT (added to the path only in Octave).
%   redbKIT calls verLessThan('matlab', ...) to choose between old/new
%   MATLAB APIs; Octave supports the "new" ones, so always return false.
tf = false;
end
