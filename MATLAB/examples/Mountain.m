classdef Mountain < grid3D
    methods(Access = public)
        function obj = Mountain
            obj = obj@grid3D(); % Call grid3D constructor
            obj.bar(0,10,0,10,0,1,0.25);
            z = sin(pi/20*obj.x).*sin(3*pi/10*obj.y).^2+0.3*rand(1,obj.nPoints);
            z = z+sin(pi/10*obj.x).*sin(3*pi/10*obj.y).^2+0.3*rand(1,obj.nPoints);
            z = z+0.1*(obj.y+obj.x)+0.1*rand(1,obj.nPoints);
            z = z + 1-0.025*((obj.x-7).^2+(obj.y-7).^2);
            obj.p(3,:) = obj.p(3,:).*(0.75*z+1);
            mdpts = obj.midpts;
            x = mdpts(1,:);
            z = mdpts(3,:);
            obj.t(5,z>=cos(pi/10*x)-0.3*x+1&z<=cos(pi/10*x)+1.5) = 2;
        end
    end
end
