classdef Newton < handle
    % A class of nonlinear solver using Newton's method
    % i.e. solves F(y) = 0
    properties
        HF; 
        % the function handle to F(y,varargin). Notice F may have other
        % input parameters
        HDF;
        % the Jacobian of the input function. Notice this is an optional
        % input parameter from the class constructor
        
    end
end