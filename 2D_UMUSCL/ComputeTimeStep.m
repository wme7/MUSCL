function dt = ComputeTimeStep(vol,wsn)
%********************************************************************************
%* This subroutine computes the explicit time-step: the minimum dt over nodes.
%*
%* ------------------------------------------------------------------------------
%*  Input: node(i).vol = Dual volume
%*         node(i).wsn = Sum of the max wave speed multiplied by the face length
%*
%* Output:         dt  = global time step
%*         node(:).dt  =  local time step
%* ------------------------------------------------------------------------------
%*
%* NOTE: Local time step is computed and stored at every node, but not used.
%*       For steady problems, it can be used to accelerate the convergence.
%*
%********************************************************************************
% Minimum time step requested by user:
dt_min = 1.0E+05;

% Local time step: dt = volume/sum(0.5*max_wave_speed*face_area).
nodes_local_dt = vol ./ wsn;

% Global time step: min(local_dt)
dt = min([dt_min;nodes_local_dt]);

end