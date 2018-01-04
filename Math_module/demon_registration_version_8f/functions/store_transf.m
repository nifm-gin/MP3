function stop=store_transf(x,optimValues,state)
% This function STORE_TRANSF is needed by the demon energy optimizer.
% It stores the Transformation field of an itteration in a global variable, needed to calculate
% transformation update energy in the next itteration.
global last_transformation_field updatemovie;
switch state
    case 'iter'
        last_transformation_field=x;
    case 'interrupt'
        clear('last_transformation_field');
    case 'init'
        last_transformation_field=x;
    case 'done'
        clear('last_transformation_field');
otherwise
end
%updatemovie=true;
updatemovie=false;
stop=false;

