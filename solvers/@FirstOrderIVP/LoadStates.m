function LoadStates(obj, u, hT)
%LoadStates load the states previously computed and advance to the latest
%state by changing CurrentIndex. The IVP object must be unsolved and
%contain matching initial condition
% WARNING! THE DERIVITAVE FUNCTION MAY NOT BE THE SAME
%   u: the states
%   t: the vector of time corresponding to each state

    % making sure it's loading into an empty instance with the same IC
    assert(obj.CurrentIndex == 1);
    assert(isequal(IC,u(1,:)));
    obj.state = u;
    obj.hT = hT;
    obj.T = hT(end);
    obj.CurrentIndex = length(t);
end

