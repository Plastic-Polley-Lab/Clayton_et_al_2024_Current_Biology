function cm = confusionchart_ke(varargin)
%CONFUSIONCHART Plot a confusion matrix
%
%   cm = CONFUSIONCHART(m) plots a confusion matrix and returns a
%   ConfusionMatrixChart object. The matrix m must be a valid confusion
%   matrix, that is, m must be a square matrix and its elements must be
%   positive integers. The element m(i,j) is the number of times an
%   observation in the ith true class was predicted to be in the jth class.
%   Each colored cell of the chart corresponds to a single element of the
%   confusion matrix.
%
%   cm = CONFUSIONCHART(m, labels) plots a confusion matrix using the
%   matrix m and the array of class names labels. Each element of labels
%   corresponds to the name of a single class. The labels array must be the
%   same size as the confusion matrix m, and can be a categorical, numeric,
%   logical, character or string vector, or a cell array of character
%   vectors.
%
%   cm = CONFUSIONCHART(trueLabels, predictedLabels) calculates and plots a
%   confusion matrix from the vectors of class labels trueLabels and
%   predictedLabels. Each element of the vectors corresponds to a single
%   observation. The vectors of class labels can be categorical, numeric,
%   logical, character or string arrays, or cell arrays of character
%   vectors.
%
%   cm = CONFUSIONCHART(___, Name, Value) specifies additional options
%   for the plot using one or more name-value pair arguments. Specify the
%   options after all other input arguments.
%
%   cm = CONFUSIONCHART(parent, ___) creates the plot in the figure,
%   panel, or tab specified by parent.

% Copyright 2017-2018 The MathWorks, Inc.

narginchk(1, Inf);

try
    % Check input args are valid, throw error here if they're not.
    [parent, model, cellOfKeyValuePairs] = iParseInput(varargin{:});
    
    % If the user hasn't specified an 'OuterPosition' value, we try and
    % replace the existing axes (if one exists). If they have specified a
    % 'Position', 'InnerPosition' or 'OuterPosition' value, we make a new
    % chart and leave any existing axes alone.
    if ~iHasPositionArg(cellOfKeyValuePairs)
        constructorFcn = @(varargin)(mlearnlib.graphics.chart.ConfusionMatrixChart(...
            varargin{:}, 'Model', model, cellOfKeyValuePairs{:}));
        
        % If the user hasn't defined a parent, parent will be empty. This
        % will be handled correctly by prepareCoordinateSystem.
        cm = matlab.graphics.internal.prepareCoordinateSystem('mlearnlib.graphics.chart.ConfusionMatrixChart', parent, constructorFcn);
    else
        % If the user hasn't defined a parent, we need to get one now.
        if isempty(parent)
           parent = gcf(); 
        end

        cm = mlearnlib.graphics.chart.ConfusionMatrixChart('Parent', parent, 'Model', model, cellOfKeyValuePairs{:});
    end
    
    fig = ancestor(cm, 'Figure');
    if isscalar(fig)
        fig.CurrentAxes = cm;
    end
catch e
    throw(e);
end

end

function [parent, model, cellOfKeyValuePairs] = iParseInput(varargin)
% Parse input to standard form. We return the parent separately (even
% though it's included in the args) so we can use it in
% matlab.graphics.internal.prepareCoordinateSystem.

parentValidator = mlearnlib.internal.confusionmatrixchart.input.ParentValidator();
syntaxInterpreter = mlearnlib.internal.confusionmatrixchart.factories.SyntaxInterpreter();

parser = mlearnlib.internal.confusionmatrixchart.input.InputParser(parentValidator, syntaxInterpreter);

[parent, model, cellOfKeyValuePairs] = parser.parse(varargin{:});
end

function tf = iHasPositionArg(cellOfKeyValuePairs)
% Returns true if the name-value pairs contain a 'Position' name of some
% sort.

propNames = cellOfKeyValuePairs(1:2:end);
tf = any(strcmpi(propNames, 'Position')) || ...
    any(strcmpi(propNames, 'InnerPosition')) || ...
    any(strcmpi(propNames, 'OuterPosition'));
end