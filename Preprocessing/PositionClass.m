classdef PositionClass < handle
    properties
        x
        y
        a
        b
        r
        alpha
    end
    
    methods
        % class constructor
        function obj = PositionClass(x,y,a,b,r,alpha)
            obj.x = x;
            obj.y = y;
            obj.a = a;
            obj.b = b;
            obj.r = r;
            obj.alpha = alpha;
        end
    end
end

    