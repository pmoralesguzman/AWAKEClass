%---------------------------------------------------------------------
% Shortcut for CreatePathClass
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
%
% P. I. Morales Guzman
% Last update: 07/01/2021
%---------------------------------------------------------------------

function cpc(varargin)

if nargin == 0
    CreatePathClass()
else
    CreatePathClass(1)
end

set(0,'DefaultAxesColorOrder',linspecer(8))
set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultFigureColormap',linspecer);
end

