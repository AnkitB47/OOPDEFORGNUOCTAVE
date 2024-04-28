classdef grid3Dpr < gridd
    %grid3Dpr
    % Creates prism element meshes.
    % (c) 2013 Uwe Prüfert
    properties(SetAccess = protected)
        basis       % Store the 2D grid object
        nPointsInElements = 6;
    end

    properties(Constant = true)
        spaceDimension = 3;
    end


    methods(Access = public)
        function obj = grid3Dpr()
            % constructor method
            % only empty call allowed
            % obj = grid3Dpr()
            switch nargin
                case 0
                    % empty object
                otherwise
                    obj.wrongNumberInputs.throwAsCaller();
            end
        end

        function obj2 = copy(obj1)
            %%
            % * copy IN:self OUT:gri3Dpr
            %
            % Hard copy method.
            obj2 = copy@gridd(obj1);
            obj2.basis = obj1.basis.copy;
        end

        % some geometries
        function unitCube(obj,h)
            switch nargin
                case 1
                    h = 0.1;
                case 2
                    %
                otherwise
            end
            obj.basis = grid2D();
            obj.basis.unitSquare(h);
            obj.extrude(obj.basis, linspace(0,1,ceil(2/h)+1));

            % Correct the Boundary numbering...
            % because the unitcube has its own system
            for k = 1:obj.nEdges
                A = obj.p(:,obj.e(1:3,k));
                if A(1,:) == zeros(1,3)
                    obj.e(5,k) = 2;
                elseif A(1,:) == ones(1,3)
                    obj.e(5,k) = 4;
                elseif A(2,:) == zeros(1,3)
                    obj.e(5,k) = 3;
                elseif A(2,:) == ones(1,3)
                    obj.e(5,k) = 5;
                elseif A(3,:) == zeros(1,3)
                    obj.e(5,k) = 1;
                elseif A(3,:) == ones(1,3)
                    obj.e(5,k) = 6;
                end
            end

        end

        function bar(obj,varargin)
            % method that meshes a 3D bar.
            % obj.bar(xmin,xmax,ymin,ymax,zmin,zmax)
            % obj.bar(xmin,xmax,ymin,ymax,zmin,zmax,hmax)
            % obj.bar(xvec,yvec,zvec)
            % xvec, ybvec zvec must be vectors with length >2.
            % This call is intentionally for using different meshwides in
            % x-y and z direction.
            % (c) 2013 by Uwe Prüfert

            obj.basis = grid2D();
            switch nargin
                case {4}
                    if any([...
                            length(varargin{1})<=2 ...
                            length(varargin{2})<=2 ...
                            length(varargin{3})<=2])
                       ME = MException('GRID3DPR:WRONGUSE',...
                            ['In the 3 argument call the arguments must',...
                            ' be vectors defining',...
                            ' the coordinates of the points of the bar.']);
                       ME.throwAsCaller();
                    else
                        % compute h for meshing the base square
                        h = min([varargin{1}(2:end)-varargin{1}(1:end-1),...
                            varargin{2}(2:end)-varargin{2}(1:end-1)]);
                    end

                    obj.basis.square(varargin{1}(1),varargin{1}(end),...
                        varargin{2}(1),varargin{2}(end),h);
                    obj.extrude(obj.basis,varargin{3});
                case {7,8}
                    switch nargin
                        case 8
                            h = varargin{7};
                        case 7
                            h = 0.1;
                        otherwise
                            obj.wrongNumberInputs.throwAsCaller();
                    end

                    obj.basis.square(varargin{1},varargin{2},varargin{3},varargin{4},h);
                    obj.extrude(obj.basis,linspace(varargin{5},varargin{6},...
                        ceil((varargin{6}-varargin{5})/h)+1));
                otherwise
                    ME = MException('GRID3DPR:WRONGNUMBEROFARGUMENTS',...
                         [' The number of arguments must be 3, 6 or 7.',...
                          'See help grid3Dpr.bar.']);
                    ME.throwAsCaller();
            end



        end

        function cylinder(obj,r,l,rmax,lmax)
            switch nargin
                case 1
                    r = 1;
                    l = 1;
                    rmax = 0.2;
                    lmax = 0.3;
                case 3
                    rmax = 0.2;
                    lmax = 0.3;
                case 4
                    lmax = 0.7*rmax;
                case 5
                    % full input list...
                otherwise
                    obj.wrongNumberInputs.throwAsCaller();
            end
            obj.basis = Circle(r,rmax);
            obj.extrude(obj.basis,linspace(0,l,ceil(l/lmax)+1));
        end

        function rail(obj,l,lmax)
            switch nargin
                case 1
                    l = 3;
                     lmax = 0.3;
                case 2
                    lmax = l/10;
                case 3
                    % full input list...
                otherwise
                    obj.wrongNumberInputs.throwAsCaller();
            end
            obj.basis = grid2D();
            obj.basis.doubleT();
            obj.basis.refineMesh();
            obj.basis.refineMesh();
            obj.extrude(obj.basis,linspace(0,l,ceil(l/lmax)+1));
        end

        function plot(obj,varargin)
            % PLOT method.
            % obj.plot([patchArguments])
            % Optinal patchArguments may all arguments that can delivered
            % to patch objects. See also help patch.
            obj.plotFaces([],'FaceColor','none',varargin{:})
            view(3);
            axis equal;
        end

        function [cval,aval,fval] = aCoefficientsMpt(obj,cc,aa,ff)
            %computes the value of the coefficients in the center of every triangle
            % 'symbolic' variables x, y  are neccesary for evaluation of string
            % objects like c = 'sin(x)' etc.
            % Restrictions
            % NOT jet implemented:
            % changeLog:
            % 03/03/2015
            % cell array input added
            % f ull matrix c input added


            % p = obj.p;
            % t = obj.t;
            np = obj.nPoints;
            nt = obj.nElements;

            switch class(cc)
%                 case 'function_handle'
%                     mpt = 1/3*(obj.p(:,obj.t(1,:))...
%                             +obj.p(:,obj.t(2,:))...
%                             +obj.p(:,obj.t(3,:)));
%                     x = mpt(1,:);
%                     y = mpt(2,:);
%                     z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
%                     try
%                         c = feval(cc,x,y,z);
%                     catch ME
%                         rethrow(ME);
%                     end
%                     cval(1,:) = c(1,1:nt));
%                     cval(2,:) = c2(ones(1,nt));
%                     cval(3,:) = c3(ones(1,nt));

                case 'double'
                    % five case
                    % (i) scalar
                    % (ii.1) 3 x 3 diagonal matrix
                    % (ii.2) 3 x 3 full matrix
                    % (iii) vector length np
                    % (iv) vector length nt
                    [rows,cols] = size(cc);
                    switch max(rows,cols)
                        case 1
                            cval = cc(ones(3,nt));
                        case 3
                           if ( norm((ones(3) - diag(3)).*cc,1) == 0)
                               % (ii.1) 3 x 3 diagonal matrix
                                c1 = cc(1,1); c2 = cc(2,2); c3 = cc(3,3);
                                cval(1,:) = c1(ones(1,nt));
                                cval(2,:) = c2(ones(1,nt));
                                cval(3,:) = c3(ones(1,nt));
                           else
                                % (ii.2) 3 x 3 full matrix
                                c1 = cc(1,1); c2 = cc(2,1); c3 = cc(3,1);
                                c4 = cc(1,2); c5 = cc(2,2); c6 = cc(3,2);
                                c7 = cc(1,3); c8 = cc(2,3); c9 = cc(3,3);
                                cval(1,:) = c1(ones(1,nt));
                                cval(2,:) = c2(ones(1,nt));
                                cval(3,:) = c3(ones(1,nt));
                                cval(4,:) = c4(ones(1,nt));
                                cval(5,:) = c5(ones(1,nt));
                                cval(6,:) = c6(ones(1,nt));
                                cval(7,:) = c7(ones(1,nt));
                                cval(8,:) = c8(ones(1,nt));
                                cval(9,:) = c9(ones(1,nt));
                           end
                        case nt
                            % the good one, nothing to do
                            if min(rows,cols)==1
                                cval = ones(3,1)*cc(:)';
                            elseif (cols==3) || (cols == 9)
                                cval = cc';
                            elseif (rows==3) || (rows == 9)
                                cval = cc;
                            else
                                obj.wrongFormat.throwAsCaller();
                            end

                        case np

                            ME = MException('grid3Dpr:aCoefficients',...
                                 ['Not jet implemented, try to use a function',...
                                 ' for the midpoints of elements']);
                            ME.throwAsCaller();

                        otherwise
                             obj.wrongFormat.throwAsCaller();
                    end
                case 'char'
                    % must be a single char symbolizing
                    % the coefficent function.
                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".
                    % Hence, the next warning can be ignored!
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
                    try
                        cval = eval(cc);
                    catch ME
                        throw(ME);
                    end
                    [rows,cols] = size(cval);
                    switch max(rows,cols)
                        case 1
                            % c is a constant like 'pi'
                            cval = cval(ones(3,nt));
                        case nt
                            cval = [cval;cval;cval];
                        otherwise
                            obj.wrongFormat.throwAsCaller();
                    end
                case 'cell'
                    % must be a 3*3 cell array which refers to the entries of a 3*3
                    % diffussion matrix. Every Entry must be double or a single char symbolizing
                    % the part of the coefficent function like in the case
                    % 'char'.

                    [rows,cols] = size(cc);
                    if (rows ~=3) || (cols ~=3)
                        obj.wrongFormat.throwAsCaller();
                    end
                    cval = zeros(9,nt);
                    cvalcell = cell(3,3);


                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".
                    % Hence, the next warning can be ignored!
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
                    try
                        for i = 1:3
                            for j = 1 : 3
                                switch class(cc{i,j})
                                    case 'char'
                                        cvalcell{i,j} = eval(cc{i,j});
                                    case 'double'
                                        cvalcell{i,j} = cc{i,j};
                                    otherwise
                                        obj.wrongClass.throwAsCaller();
                                end
                            end
                        end
                    catch ME
                        throw(ME);
                    end
                        for i = 1:3
                            for j = 1 : 3
                                linIndex = 3*(i-1)+j;
                                cvalij = cvalcell{i,j};
                                cvalij = cvalij(:);
                                [rows,cols] = size(cvalij);
                                switch max(rows,cols)
                                    case 1
                                        % c is a constant like 'pi'
                                        cval(linIndex,:) = cvalij(ones(1,nt));
                                    case nt
                                        cval(linIndex,:) = cvalij';
                                    otherwise
                                        obj.wrongFormat.throwAsCaller();
                                end
                            end
                        end
                otherwise
                    obj.wrongClass.throwAsCaller();
            end
            % repeat code from case 'c'...
            switch class(aa)
                case 'double'
                    [rows,cols] = size(aa);
                    switch max(rows,cols)
                        case 1
                            aval = aa(ones(1,nt));
                        case np
                            aval = obj.point2Center(aa);
                            aval = aval(:)';
                        case nt
                            % the good one, nothing to do
                            aval = aa(:)';
                        otherwise
                            obj.wrongSize.throwAsCaller();
                    end
                case 'char'
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
                    try
                        aval = eval(aa);
                    catch ME
                        throw(ME);
                    end
                    [rows,cols] = size(aval);
                    switch max(rows,cols)
                        case 1
                            aval = aval(ones(1,nt));
                        case nt
                            % work already done
                        otherwise
                            finiteElements.wrongFormat.throwAsCaller();
                    end
                case 'cell'
                    error('Sorry, Cell array input not jet implemented!');
                otherwise
                    obj.wrongClass.throwAsCaller();
            end
            switch class(ff)
                case 'double'
                    [rows,cols] = size(ff);
                    switch max(rows,cols)
                        case 1
                            fval = ff(ones(1,nt));
                        case np
                            fval =  obj.point2Center(ff);
                        case nt
                            % the good one, nothing to do
                            fval = ff(:)';
                        otherwise
                            obj.wrongSize.throwAsCaller();
                    end
                case 'function_handle'
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));
                    try
                        fval = feval(ff,x,y,z);
                    catch ME
                        throw(ME);
                    end
                case 'char'
                    mpt = 1/3*(obj.p(:,obj.t(1,:))...
                            +obj.p(:,obj.t(2,:))...
                            +obj.p(:,obj.t(3,:)));
                    x = mpt(1,:);
                    y = mpt(2,:);
                    z =  0.5*(obj.p(3,obj.t(1,:))+ obj.p(3,obj.t(4,:)));

                   try
                        fval = eval(ff);
                    catch ME
                        throw(ME);
                    end

                    [rows,cols] = size(fval);
                    switch max(rows,cols)
                        case 1
                           fval = fval(ones(1,nt));
                        case nt
                            % work already done
                        otherwise
                            obj.wrongFormat.throwAsCaller();
                    end
                case 'cell'
                    error('Sorry, Cell array input makes here no sense.');
                otherwise
                    obj.wrongClass.throwAsCaller();
            end
        end

        function [bvalvec] = cCoefficients(g,b)
            %computes the value of the coefficient b in the center of every triangle
            % coefficent can be a vector of dim 2 x 1, or a cell array of dim 2 x 1

            % 'symbolic' variables x, y  and t are neccesary for evaluation of string
            % objects like c = 'sin(x)' etc.
            p = g.p;
            t = g.t;

            x = p(1,:);
            y = p(2,:);
            z = p(3,:);

            % number of points and triangles
            n = g.nPoints;
            nt = g.nElements;

            % check the class and  size of b
            if ~((max(size(b))>=2)&&(min(size(b))>=1))
                ME = MException('ccoefficients:wrongCoefficientDefinition',...
                    ' b must be a vector');
                ME.throwAsCaller();
            end
            % two cases:
            % cell-array - entries are strings, or doubles
            if isa(b,'cell')
                b1 = b{1};
                b2 = b{2};
                b3 = b{3};
            elseif isa(b,'double')
                b1 = b(1,:);
                b2 = b(2,:);
                b3 = b(3,:);
            else
                ME = MException('ccoefficients:wrongCoefficientDefinition',...
                    ' Wrong coefficient definition');
                ME.throwAsCaller();
            end


            if isa(b1,'function_handle') || isa(b1,'inline')
                bval = feval(b1,x,y);
            elseif isa(b1,'char'),
                bval = eval(b1).*ones(1,n);
            elseif isa(b1,'numeric')
                if length(b1)==n
                    % c vektor and defined in p
                    bval = b1;
                elseif length(b1)==1,
                    % skalar
                    bval = b1*ones(nt,1);
                elseif length(b1)==nt,
                    bval = b1;
                else
                    ME = MException('ccoefficients:wrongSize','wrong sized b(1)');
                    ME.throwAsCaller();
                end
            elseif isa(b1,'inline')
                 bval = b(x,y);
            else
                ME = MException('ccoefficients:wrongSize','wrong formated b(1)');
                ME.throwAsCaller();
            end
            if length(bval) == nt
                % b(1) is a vektor and defined in center of mass of triagle
            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = g.point2Center(bval);
            end
            dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
            % first column of bvalvec
            bvalvec = bval;

            if isa(b2,'function_handle') || isa(b2,'inline')
                bval = feval(b2,x,y);
            elseif isa(b2,'char'),
                bval = eval(b2).*ones(1,n);
            elseif isa(b2,'numeric')
                if length(b2)==n
                    % c vektor and defined in p
                    bval = b2;
                elseif length(b2)==1,
                    % skalar
                    bval = b2*ones(nt,1);
                elseif length(b2)==nt,
                    bval = b2;
                else
                    ME = MException('ccoefficients:wrongSize','wrong sized b(1)');
                    ME.throwAsCaller();
                end
            elseif isa(b2,'inline')
                 bval = b(x,y);
            else
                ME = MException('ccoefficients:wrongSize','wrong formated b(1)');
                ME.throwAsCaller();
            end

            if length(bval) == nt
                % b(1) is a vektor and defined in center of mass of triagle

            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = g.point2Center(bval);
            end
            dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
            bvalvec = [bvalvec,bval];

            if isa(b3,'function_handle') || isa(b2,'inline')
                bval = feval(b2,x,y);
            elseif isa(b3,'char'),
                bval = eval(b2).*ones(1,n);
            elseif isa(b3,'numeric')
                if length(b3)==n
                    % c vektor and defined in p
                    bval = b3;
                elseif length(b3)==1,
                    % skalar
                    bval = b3*ones(nt,1);
                elseif length(b3)==nt,
                    bval = b3;
                else
                    ME = MException('ccoefficients:wrongSize','wrong sized b(1)');
                    ME.throwAsCaller();
                end
            elseif isa(b3,'inline')
                 bval = b(x,y);
            else
                ME = MException('ccoefficients:wrongSize','wrong formated b(1)');
                tME.throwAsCaller();
            end
            if length(bval) == nt
                % b(1) is a vektor and defined in center of mass of triagle

            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = g.point2Center(bval);
            end

            dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
            bvalvec = [bvalvec,bval]';

        end

        function b = isPointInDomain(obj,pt)
            % isPointInTetraeder -  gives back a bool vector of length nElements.
            % b ist true if the point pt is i the element of a given
            % mesh object.
            % b = isPointInElement(obj,pt)
            % Algortihm: Transform the point wrt the transformation of
            % the elements into the unit element  and
            % decide on basis of the relation pt_x < 0 pt_y < 0 ...
            p = obj.p;
            t = obj.t;

            p1 = p(:,(t(1,:)));
            p2 = p(:,(t(2,:)));
            p3 = p(:,(t(3,:)));
            p4 = p(:,(t(4,:)));

            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:);
            x41 = p4(1,:)-p1(1,:);

            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:);
            y41 = p4(2,:)-p1(2,:);

            z21 = p2(3,:)-p1(3,:);
            z31 = p3(3,:)-p1(3,:);
            z41 = p4(3,:)-p1(3,:);


            J =  ( (x21.*y31-x31.*y21).*z41...
                -(x21.*y41-x41.*y21).*z31...
                +(x31.*y41-x41.*y31).*z21);

            % transform the point into the unit element and
            % decide

            % to compute     the inverse of linear transformation
            xi_x = (y31.*z41-y41.*z31)./J;
            eta_x = (y41.*z21-y21.*z41)./J;
            zeta_x =  (y21.*z31-y31.*z21)./J; % zero for prisms

            xi_y = (x41.*z31-x31.*z41)./J ;
            eta_y = (x21.*z41-x41.*z21)./J;
            zeta_y = (x31.*z21-x21.*z31)./J; % zero for prisms

            xi_z = (x31.*y41-x41.*y31)./J;
            eta_z = (x41.*y21-x21.*y41)./J;
            zeta_z = (x21.*y31-x31.*y21)./J;

            b = isPointInAnyTriangle3DPR(pt,p1,xi_x,...
                                         eta_x,xi_y,...
                                         eta_y,xi_z,...
                                         eta_z,zeta_z) ;
            % If we use the MEX version, we must convert it into logical.
            b = logical(b);
%
            % There is a MEX function isPointInAnyTriangle
            % If you wish to use it, and you have compiled it for your OS,
            % comment out the following function. Note that the speed up is
            % not dramatically.

            function b = isPointInAnyTriangle3DPR(pt,p1,xi_x,...
                                                  eta_x,xi_y,...
                                                  eta_y,xi_z,...
                                                  eta_z,zeta_z)

                nP = size(pt,2);
                b = false(1,nP);

                for k = 1:nP
                    ptx = pt(1,k)-p1(1,:);
                    pty = pt(2,k)-p1(2,:);
                    ptz = pt(3,k)-p1(3,:);

                    ptux = xi_x.*ptx+xi_y.*pty+xi_z.*ptz;
                    ptuy = eta_x.*ptx+eta_y.*pty+eta_z.*ptz;
                    ptuz = zeta_z.*ptz;

                    b(k) = any((sum([ptux;ptuy])<=1 &  abs(ptuz)<=1 &  ptux>=0 &  ptuy>=0 ));
                end
            end
        end

        function [sidelength,area] = sideLengthAndArea(obj)
            % abstract defined in gridd
            % Must be implemented!
        end

        function extrude(obj,grid2d,z)
        %extrude(obj,grid2d,z)
            %  Extrudes a 2D geometry in z direction. grid2d must be a valid GRIDD2D
            %  object defining a 2D geometry object. The result is a
            %  6-NODE WEDGE ELEMENT (aka prism element) mesh.

            switch class(z)
                case 'grid1D'
                    z = unique(z.p);
                case 'double'
                    z = unique(z);
                otherwise
                    ME = MException(...
                         'GRID3DPR:wrongClass',...
                         ['The second argument was a',...
                         class(z), ' Allowed is drid1D or double']);
                    ME.throwAsCaller();
            end

                % To be really sure and make it fools save ;-)
            if ~isa(grid2d,'grid2D')|| min(size(grid2d.p)) ~= 2
                obj.wrongClass.throwAsCaller();
            else
                obj.basis = grid2d;
            end
            nSlices = length(z);
            nPts = obj.basis.nPoints;
            nTri =  obj.basis.nElements;
            if num2str(nTri*(nSlices-1))>1e6
                fprintf(['Warning: A mesh with ',num2str(nTri*(nSlices-1)),...
                    ' elements and ',num2str(nSlices*nPts),...
                    ' points will now be generated.\n']);
            end
            obj.p = zeros(3,nPts*nSlices);
            block = expandFirstSlice(obj.basis.t,nPts);
            obj.t = block;

            % the edges
            % initialize with triangel data from 2D object
            obj.e = obj.basis.t(1:3,:);
            % base gets the subdomain numbers from 2D object
            obj.e(5,1:nTri) = obj.basis.t(4,:) ;
            nBdConds = max(obj.e(5,:));
            EdgeBlock = expandBoundary(obj.basis.e,nPts,nBdConds);
            obj.e = [obj.e EdgeBlock];

            % wedge and points-loop
            for k = 1:nSlices-2
                block(1:6,:) = block(1:6,:)+nPts; % !!!!! Only point index, not domain index
                EdgeBlock(1:4,:) = EdgeBlock(1:4,:)+nPts;
                obj.t = [obj.t block];
                obj.p(:,nPts*(k-1)+1:nPts*k) = ...
                    [obj.basis.p;z(k)*ones(1,nPts)];
                obj.e = [obj.e [EdgeBlock(1:4,:);EdgeBlock(5,:)]];
            end
            % points-loop two more
            for k = nSlices-1:nSlices
                % build upo the points
                obj.p(:,nPts*(k-1)+1:nPts*k) = ...
                    [obj.basis.p;z(k)*ones(1,nPts)];

            end
            obj.e =  [obj.e [obj.basis.t(1:3,:)+nPts*(nSlices-1);...
               zeros(1,size(obj.basis.t,2));...
               obj.basis.t(4,:)+max(obj.e(5,:))]];


            % local block definition function, different from tetrahedra
            function block = expandFirstSlice(t,nPts)
                block = zeros(7,size(t,2));
                for kb = 1:size(t,2)
                    block(:,kb) =  [t(1,kb);...
                        t(2,kb);...
                        t(3,kb);...
                        nPts+[t(1,kb);...
                                t(2,kb);...
                                t(3,kb)];...
                        t(4,kb)]; % !!!! Add number of triangel from 2D grid
                end
            end

            function block = expandBoundary(e,nPts,nBdConds)
                % expands the boundary conditions to the mantle
                % of the cylinder
                    block  = [e(1:2,:);e(2:-1:1,:)+nPts];
                    block = [block;e(5,:)+nBdConds];

            end
        end

        function refineMesh(obj,varargin)
            % obj.refineMesh() uniform refinement
            % obj.refineMesh(indxXY,indxZ)
            % Refines all elements in indxXY in the X-Y- plane and all
            % layeres wrt. the Z axis in indxZ.
            z = unique(obj.p(3,:));
            switch nargin
                case 1
                    z = unique([z,0.5*(z(1:end-1)+z(2:end))]);
                    obj.basis.refineMesh;
                    obj.extrude(obj.basis,z);
                case 3
                    obj.basis.refineMesh(varargin{1});
                    zr = 0.5*(z(1:end-1)+z(2:end));
                    zr = zr(varargin{2});
                    z = unique([z,zr]);
                    obj.extrude(obj.basis,z);
                otherwise
            end

        end

        function rotateMesh(obj,varargin)
            switch length(varargin)
                case 4
                    rotateMesh@gridd(obj,varargin{:})
                    if varargin{1}~=0 ||varargin{2}~=0
                        warning('You try to rotate a prism grid around x or y axis.')
                    end
                case 7
                    rotateMesh@gridd(obj,varargin{:})
                    if varargin{5}~=0 ||varargin{6}~=0
                        warning('You try to rotate a prism grid around x or y axis.')
                    end
            end
        end
    end

    methods(Static)
        function[qval,gval,hval,rval] = bCoefficients(p,b)
            % compute the boundary coefficients
            %
            % q,g,h,r are srings defining everything that can be
            % evaluated by eval. The
            % independent variables must be named by x, y, z
            %
            % are three or four points
            x = sum(p(1,:))/size(p,2); %#ok<*NASGU>
            y = sum(p(2,:))/size(p,2);
            z = sum(p(3,:))/size(p,2);

            m = b(2);
            qval = 0;
            gval = 0;
            hval = 0;
            rval = 0;
            lengthq = b(3);
            lengthg = b(4);
            if m == 0 % only Neumann BCs
                qval = eval(char(b(5:5+lengthq-1)));
                gval = eval(char(b(5+lengthq:5+lengthq+lengthg-1)));
            else % only Dirichlet BCs
                lengthh = b(5);
                lengthr = b(6);
                char(b(9:9+lengthh-1));
                hval = eval(char(b(9:9+lengthh-1)));
                rval = eval(char(b(9+lengthh:9+lengthh+lengthr-1)));
            end
        end
    end

    methods(Access = public)
        function plotFaces(obj,f,varargin)
            %%
            % * plotFaces In:self,[double,[char,char]]
            %
            % Plots the 3D object. Thi first optionla argument is a vector
            % of length nPoints or number of nEdges. The surface of the
            % object will then colored depending on this data. Following
            % arguments are optional and control the plot.
            %
            % Call
            % obj.plotFaces
            % obj.plotFaces(y)
            % obj.plotFaces(y,[args]))
            %
            % plot grid faces, or plot y over faces.
            % triangle

            switch nargin
                case 1
                    f = [];
                otherwise
                    if length(f) == obj.nEdges
                        % data on barycenters of faces
                        f = repmat(f(:)',4,1);
                    end
            end

            %triangle elements
            indx = find(obj.e(4,:)==0);
            pts = (obj.e(1:3,indx)); % triangle points
            if isempty(f)
                % Fill with gray color.
                c = [  0.9   0.9    0.8];
            elseif isvector(f)
                c = f(pts);
            else
                c = f(1:3,indx);
            end

            x = reshape(obj.p(1,obj.e(1:3,indx)),3,length(indx));
            y = reshape(obj.p(2,obj.e(1:3,indx)),3,length(indx));
            z = reshape(obj.p(3,obj.e(1:3,indx)),3,length(indx));

            patch(x,y,z,c,varargin{:});

            % square elements
            indx = find(obj.e(4,:)>0);
            pts = (obj.e(1:4,indx));
            if isempty(f)
                %
            elseif isvector(f)
                c = f(pts);
            else
                c = f(1:4,indx);
            end

            x = reshape(obj.p(1,obj.e(1:4,indx)),4,length(indx));
            y = reshape(obj.p(2,obj.e(1:4,indx)),4,length(indx));
            z = reshape(obj.p(3,obj.e(1:4,indx)),4,length(indx));

            patch(x,y,z,c,varargin{:});

            view(3)
            axis equal
            if nargin>2 && ~isempty(f)
                colorbar
            end
        end

        function plotIso(obj,u,lev,varargin)
        % iso-surfaceplot
        % obj.plotIso(y,iso,options)
        % Plots the iso-surface of the data y to the level iso
        % Optional arguments are NGP (Number of Grid Points of the iso
        % surface) and the FontSize of the plot.
        % This method bases on the request and the prototype code
        % of Hannes Uecker. Thanks!
        % (C) 2014 by Uwe Prüfert


            ng = 20;
            fs = 14;
            color = [0.6 0.3 0.8
                     0.8 0.8 0.8
                     0.5 0.5 1
                     0.6 0.6 1
                     0.2 0.7 0.8];

            if all(lev>max(u)) && all(lev < min(u))
                fprintf('Empty level set\n');
            end

            gp=obj.p;

            x1=min(gp(1,:));
            x2=max(gp(1,:));
            y1=min(gp(2,:));
            y2=max(gp(2,:));
            z1=min(gp(3,:));
            z2=max(gp(3,:));

            xv=linspace(x1,x2,ng);
            yv=linspace(y1,y2,ng);
            zv=linspace(z1,z2,ng);

            [X,Y,Z]=meshgrid(xv,yv,zv);

            up=obj.p3interpol(X,Y,Z,u,gp(1,:),gp(2,:),gp(3,:));

            if length(lev)>5
                warning('PLOTUTIL3D:TOOMANYLEVELS',...
                    'The number of iso surfaces to plot is restricted to five.')
            end
            hold on
            for k=1:min(length(lev),5)
                ip=patch(isosurface(X,Y,Z,up,lev(k)));
                isonormals(X,Y,Z,up,ip);
                set(ip,'FaceColor',color(k,:),'EdgeColor','none');
            end

            view(3); axis([x1 x2 y1 y2 z1 z2]);
            camlight;
            lighting phong;
            xlabel('x','FontSize',fs);
            ylabel('y','FontSize',fs);
            zlabel('z','FontSize',fs);
            grid on;
            title(['Level = ' mat2str(lev,3)],'FontSize',fs);
            set(gca,'FontSize',fs);
            plot3(x1,y1,z1);
            plot3(x2,y2,z2);
            axis equal
            hold off
            if nargin>3
                % additional options
                set(ip,varargin{:});
            end
        end

        function cutawayPlot(obj,varargin)
        % cutawayPlot - cuts away all points x<X, y<Y z<Z and
        % plots the object.
        % obj.cutawayPlot(X,Y,Z)
        % Two of the arguments X,Y,Z may be empty.
        % obj.cutawayPlot(X,Y,Z,y) plots the object with the data y
        % (c) 2014 by Uwe Prüfert


            % make a copy of all data
            p = obj.p;

            % we do not need the edge number information here...
            e = obj.e(1:4,:);
            t = obj.t;

            % find all point with x<=X etc...

            switch length(varargin)
                case 3
                    X = varargin{1};
                    Y = varargin{2};
                    Z = varargin{3};
                    X2 = [];
                    Y2 = [];
                    Z2 = [];
                    y = [];
                case 4
                    X = varargin{1};
                    Y = varargin{2};
                    Z = varargin{3};
                    X2 = [];
                    Y2 = [];
                    Z2 = [];
                    y = varargin{4};
                case 6
                    X = varargin{1};
                    Y = varargin{2};
                    Z = varargin{3};
                    X2 = varargin{4};
                    Y2 = varargin{5};
                    Z2 = varargin{6};
                case 7
                   X = varargin{1};
                   Y = varargin{2};
                   Z = varargin{3};
                   X2 = varargin{4};
                   Y2 = varargin{5};
                   Z2 = varargin{6};
                   y = varargin{7};
                otherwise
                    obj.wrongNumberInputs.throwAsCaller();
            end
            if ~isempty(X)
                indX = find(obj.p(1,:)<X);
            else
                indX = [];
            end
            if ~isempty(Y)
                indY = find(obj.p(2,:)<Y);
            else
                indY = [];
            end
            if ~isempty(Z)
                indZ = find(obj.p(3,:)<Z);
            else
                indZ = [];
            end
            if ~isempty(X2)
                indX2 = find(obj.p(1,:)>X2);
            else
                indX2 = [];
            end
            if ~isempty(Y2)
                indY2 = find(obj.p(2,:)>Y2);
            else
                indY2 = [];
            end
            if ~isempty(Z2)
                indZ2 = find(obj.p(3,:)>Z2);
            else
                indZ2 = [];
            end
            indX = unique([indX indY indZ indX2 indY2 indZ2]);
            if isempty(indX)
                ME = MException('grid3Dpr:cutawayPlot:EmptySet',...
                     'The cut is not well defined, check the arguments.');
                ME.throwAsCaller();

            end

            % remove points
            numberOfRemovedPoints = length(indX);
            fprintf([num2str( numberOfRemovedPoints),' points are in the cutaway.\n'])

            for k = indX
                % look for removed points in elements and set the
                % position of removed points to zero.
                % indX(k) found in r-th row and c-th column

                [r,c] = find(e==k);

                for k2 = 1:length(r)
                    e(r(k2),c(k2)) = 0;
                end
                [r,c] = find(t==k);

                for k2 = 1:length(r)
                    t(r(k2),c(k2)) = 0;
                end
            end

            % now all all zeros t can be removed and
            % all t with less than 3 non-zeros
            % all three and four non-zeros become new edges
            k = 1;
            % remove deformed edges
            while true
                if k>size(e,2)
                    break;
                end
                switch length(find(e(:,k)>0))
                    case {3,4}
                        % edge element triangle or
                        % edge element rectangle
                        % stays in e so we skip the removing
                        k = k+1;
                    otherwise % 0 1 2
                        % to removed: patches with removed points
                        % do not increase k because k + 1 "jumps" to k
                        % "automaticly"
                        e(:,k) = [];
                end
            end

            % add former elements (now with three or four points insted
            % of six to edges
            for k = 1:obj.nElements
                switch length(find(t(:,k)>0))
                    case 3
                        % edge element triangle or
                        % edge element rectangle
                        % stays in e so we skip the removing
                        e = [e ,[t((t(:,k)>0),k);0]];
                    case 4
                        % edge element triangle or
                        % edge element rectangle
                        % stays in e so we skip the removing
                        id = find(t(:,k)>0);
                        % correct the order of points in the "cut in x or
                        % cut in y" case
                        e = [e ,t(id([1 2 4 3]),k)];
                    otherwise % 0 1 2
                        %  do nothing
                end
            end

            % create a grid3Dpr object and fill it with the OLD points and
            % NEW edges. The t is a dummy, becaus grid3Dpr.display needs a
            % nonempty t. It 's only to prevent "mistyrious errors"
            p_save = obj.p;
            e_save = obj.e;
            t_save = obj.t;


            obj.p = p;
            obj.e = e;
            obj.t = 1;


            % use the plotFaces
            if isempty(y)
                obj.plotFaces();
            else
                obj.plotFaces(y);
            end
            obj.p = p_save;
            obj.e = e_save;
            obj.t = t_save;
        end

        function plotSlices(obj,u,xslice,yslice,zslice,varargin)
        % plotSlices - plots slices of a 3D function over 3D grid object.
        % plotSlices(obj,u)
        % plotSlices(obj,u,xslice,yslice,zslice)
        % plotSlices(obj,u,xslice,yslice,zslice[,options])
        % xslice,yslice,zslice are vectors given the coordinates of
        % the slices.
        % Every of the arguments xslice,yslice,zslice can be empty.
        % In this case it will be set to a mean value.
        % This method bases on the request and the prototype code
        % of Hannes Uecker. Thanks!
        % (C) 2014 by Uwe Prüfert

            x1 = min(obj.p(1,:)); x2 = max(obj.p(1,:));
            y1 = min(obj.p(2,:)); y2 = max(obj.p(2,:));
            z1 = min(obj.p(3,:)); z2 = max(obj.p(3,:));

            %
            ng = 40;
            fs = 14;

            switch nargin
                case {1,2}
                    %only u given.
                    xslice = 0.5*(x1+x2);
                    yslice = 0.5*(y1+y2);
                    zslice = 0.5*(z1+z2);
                otherwise
                    % full set of parameters given
            end

            if isempty(xslice)&&isempty(yslice)&&isempty(zslice)
                xslice = 0.5*(x1+x2);
                yslice = 0.5*(y1+y2);
                zslice = 0.5*(z1+z2);
            end
            xv=linspace(x1,x2,ng); yv=linspace(y1,y2,ng); zv=linspace(z1,z2,ng);

            xv = unique([xv , xslice]);
            yv = unique([yv , yslice]);
            zv = unique([zv , zslice]);

            [X,Y,Z] = meshgrid(xv,yv,zv);
            up = obj.p3interpol(X,Y,Z,u,obj.p(1,:),obj.p(2,:),obj.p(3,:));

            sl = slice(X,Y,Z,up,xslice,yslice,zslice);
            view(3); axis([x1 x2 y1 y2 z1 z2]);
            xlabel('x','fontsize',fs); ylabel('y','fontsize',fs);
            zlabel('z','fontsize',fs);set(gca,'FontSize',fs)
            set(sl,'EdgeColor','none','FaceColor','Interp');
            if ~isempty(varargin)
                set(sl,varargin{:});
            end
            colorbar
            hold on
            plot3(x1,y1,z1);
            plot3(x2,y2,z2);
            axis equal
        end

        function identifyBoundarySegment(obj,bdno)
            % identifyBoundary - plots the boundary segments in
            % colors.
            % obj.identifyBoundary every segment has its own color,
            % starting from blue over green to red.
            % obj.identifyBoundary(BDNR) plots the BDNR-th boundary
            % segment in yellow, all other boundary segments in blue.
            % (c) 2014 by Uwe Prüfert
            if nargin == 2
                c = double(obj.e(5,:)==bdno);
            else
                n = ceil(obj.nBoundarySegments/3);
                if n >3
                    warning({['Large number of boundary segments',...
                        ' could be unsatisfactory in this way.'];...
                        'Use identifyBoundarySegment(SegmentNumber) insted'},...
                        'Warning');

                end
                for k = 1:obj.nBoundarySegments
                     subplot(n,3,k)
                     title(['Boundary segment no ', num2str(k),])
                     obj.identifyBoundarySegment(k)
                     axis off
                end
                return
            end
            obj.plotFaces(c,'LineStyle','none');
            colorbar off;
        end

        function indx = pointToElementIndex(obj,pt)
            % pointToTriangleIndex-  gives back the index of the triangle
            % containing the point pt pt must be a vector of length
            % two.
            % indx = isPointInTriangle(obj,pt)

            b =  obj.isPointInElement(pt);
            indx = find(b);
        end
    end

    methods(Static,Access = public)
        function un = p3interpol(xn,yn,zn,u,x,y,z)
            % interpolate u def. on x,y,z to new mesh given in xn,yn,zn
            xv = reshape(x, size(x,1)*size(x,2), 1);
            yv = reshape(y, size(y,1)*size(y,2), 1);
            zv = reshape(z, size(z,1)*size(z,2), 1);
            uv = reshape(u, size(u,1)*size(u,2), 1);
            Fu = scatteredInterpolant(xv,yv,zv,uv);
            un = Fu(xn,yn,zn);
        end
    end
end

