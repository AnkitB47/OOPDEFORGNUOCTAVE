classdef RoomHeating < Elliptic
    %RoomHeating Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wood = 1e-2;
        concrete = 1e-6
        glass = 1e-1
        heater = 1e1
        temperaturBuilding = 20
        temperaturOutside = 10
    end
    
    methods(Access = public)        
        function obj = RoomHeating
            obj = obj@Elliptic;
            obj.fem = Lagrange13D;
            obj.grid = Room(0.125);
            obj.setBoundaryConditions(...
                'Robin',{obj.heater,obj.heater*25},...
                'Robin',{obj.concrete,obj.concrete*obj.temperaturBuilding},...
                'Robin',{obj.concrete,obj.concrete*obj.temperaturOutside},...
                'Robin',{obj.wood,obj.wood*obj.temperaturBuilding},...
                'Robin',{obj.glass,obj.glass*obj.temperaturOutside},...
                'Robin',{obj.glass,obj.glass*obj.temperaturBuilding});
            obj.initialize(1e-1,0,0);            
        end               
    end
end

