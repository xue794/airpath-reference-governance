function [upMode,distMode]=viD_upstream_mode(mode)
%VID_UPSTREAM_MODE Map VI-D case names to VI-A upstream and torque residual.
switch lower(mode)
    case 'aggressive'
        upMode='aggressive'; distMode='none';
    case 'aggressive_disturbed'
        upMode='aggressive'; distMode='disturbed';
    case 'wrap'
        upMode='wrap'; distMode='none';
    case 'varying'
        upMode='varying'; distMode='none';
    otherwise
        error('Unknown VI-D mode: %s',mode);
end
end
