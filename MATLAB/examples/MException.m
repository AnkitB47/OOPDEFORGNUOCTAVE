########################################################################
##
## Copyright (C) 2022 The Octave Project Developers
##
## See the file COPYRIGHT.md in the top-level directory of this
## distribution or <https://octave.org/copyright/>.
##
## This file is part of Octave.
##
## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.
##
########################################################################

classdef MException

  ## -*- texinfo -*-
  ## @deftypefn {} {@var{ME} =} MException (@var{id}, @var{template}, @dots{})
  ## Create an object of the MException class that stores error message information.
  ##
  ## The format of the input arguments matches the one from the @code{error}
  ## function.
  ## @seealso{error, throw, rethrow, throwAsCaller, getReport}
  ## @end deftypefn

  properties (GetAccess = public, SetAccess = private)
    identifier = "";
    message = "";
    stack = struct([]);
    cause = {};
    Correction = [];
    ## Universal error instances of finiteElements > abstractClasses
  end

  properties (Access = private, Hidden = true)
    hasBeenCaught = false;
  end

  methods (Access = public)

    function this = MException (id, template, varargin)

      if ((nargin == 1 || nargin == 2) && isstruct (id))
        # Matlab does not have this constructor. Used internally by
        # last() and by the evaluator in constructing exception
        # objects in CATCH blocks.
        this.identifier = id.identifier;
        this.message = id.message;
        this.stack = id.stack;
        if (nargin == 2)
          this.hasBeenCaught = true;
        end
      else
        if (nargin < 2 || nargout > 1)
          print_usage ("MException");
        end
        this.identifier = id;
        this.message = sprintf (template, varargin{:});
      end
      if (! ischar (this.identifier)
          || (! isempty (this.identifier) && ! isrow (this.identifier)))
        error ("MException: Identifier must be a row string");
      end
      if (! ischar (this.message)
          || (! isempty (this.message) && ! isrow (this.message)))
        error ("MException: Message must be a row string");
      end
      if (! isempty (this.identifier)
          && isempty (regexp (this.identifier,
                              "^([A-Za-z][A-Za-z0-9_\-]*)(:[A-Za-z][A-Za-z0-9_\-]*)+$")))
        error ("MException: Invalid identifier");
      end
    end

    function throw (this)

      ## -*- texinfo -*-
      ## @deftypefn {} {} throw (@var{ME})
      ## Issue an error from the information within an MException object.
      ##
      ## @seealso{MException, rethrow, throwAsCaller}
      ## @end deftypefn

      rethrow (as_struct (this));
    end

    function rethrow (this)

      ## -*- texinfo -*-
      ## @deftypefn {} {} rethrow (@var{ME})
      ## Reissue a previously caught exception from a try/catch statement.
      ##
      ## @seealso{MException, throw, throwAsCaller}
      ## @end deftypefn

      if (this.hasBeenCaught)
        throw (this);
      else
        error ("MException: Can only rethrow a caught exception");
      end
    end

    function throwAsCaller (this)

      ## -*- texinfo -*-
      ## @deftypefn {} {} throwAsCaller (@var{ME})
      ## Issue an error from the information within an MException object.
      ##
      ## @seealso{MException, throw, rethrow}
      ## @end deftypefn

      struct_err = as_struct (this);
      if (numel (struct_err.stack) > 0)
        struct_err.stack = struct_err.stack(2:end);
      end
      rethrow (struct_err);
    end

    function obj = addCause (this, ME)

      ## -*- texinfo -*-
      ## @deftypefn {} {@var{ME_C} =} addCause (@var{ME}, @var{C})
      ## Store another cause in an MException object.
      ##
      ## @var{C} is a scalar MException object.
      ## @seealso{MException, addCorrection}
      ## @end deftypefn

      obj = this;
      if (! isa (ME, "MException") || ! isscalar (ME))
        error ("MException: Cause must be a MException scalar object");
      end
      obj.cause{end+1} = ME;
    end

    function obj = addCorrection (this, C)

      ## -*- texinfo -*-
      ## @deftypefn {} {@var{ME_C} =} addCorrection (@var{ME}, @var{C})
      ## Store another correction in an MException object.
      ##
      ## @var{C} is a scalar object of class matlab.lang.correction.Correction.
      ## @seealso{MException, addCause}
      ## @end deftypefn

      obj = this;
      if (! isa (C, "matlab.lang.correction.Correction") || ! isscalar (C))
        error ("MException: Correction must be a Correction scalar object");
      end
      obj.Correction(end+1) = C;
    end

    function report = getReport (this, type, varargin)

      ## -*- texinfo -*-
      ## @deftypefn {} {@var{report} =} getReport (@var{ME})
      ## @deftypefnx {} {@var{report} =} getReport (@var{ME}, @var{type}, @dots{})
      ## Return the error message stored in an MException object.
      ##
      ## @var{type} is a string describing the verbosity of the error message.
      ## Default @qcode{basic} returns the error message only, while
      ## @qcode{extended}, the default, wil also return the stack of calling
      ## functions.
      ##
      ## Extra options can be specified as pairs of parameters/values. Available
      ## options are:
      ##   "hyperlinks": "on", "off", "default"
      ## @end deftypefn

      hyperlinks = "default";
      if (nargin < 2)
        type = "extended";
      end
      if (numel (varargin) > 0)
        if (mod( numel (varargin), 2))
          error ("MException: Parameter missing value");
        end
        for i = 1:2:numel (varargin)
          if (! ischar (varargin{i}))
            error ("MException: Expecting pairs of parameters/values");
          end
          switch (lower (varargin{i}))
            case "hyperlinks"
              hyperlinks = lower (varargin{i+1});
              if (! ismember (hyperlinks, {"on", "off", "default"}))
                error ("MException: hyperlinks parameter must be one of 'on', 'off', 'default'");
              end
            otherwise
              error ("MException: Unknown parameter");
          endswitch
        endfor
      end
      switch (type)
        case "basic"
          report = this.message;
        case "extended"
          report = {this.message};
          if (numel (this.stack))
            report{end+1} = "error: called from";
            for i = 1:numel (this.stack)
              report{end+1} = sprintf ("    %s at line %d column %d",  ...
                this.stack(i).name, this.stack(i).line, this.stack(i).column);
            endfor
          end
          report = [report; repmat({newline()}, 1, numel (report))];
          report = [report{1:end-1}];
        otherwise
          error ("MException: Report type must be 'basic' or 'extended'");
      endswitch
    end

    function newobj = horzcat (varargin)
      newobj = varargin{1};
      sz = cumsum (cellfun (@(x) size (x, 2), varargin));
      for i = 2:numel (varargin)
        newobj(:,sz(i-1)+1:sz(i)) = varargin{i};
      endfor
    end

    function newobj = vertcat (varargin)
      newobj = varargin{1};
      sz = cumsum (cellfun (@(x) size (x, 1), varargin));
      for i = 2:numel (varargin)
        newobj(sz(i-1)+1:sz(i),:) = varargin{i};
      endfor
    end

  end

  methods (Access=private)

    function struct_err = as_struct (this)
      if (! isscalar (this))
        error ("MException: Expect a scalar input");
      end
      if (isstruct (this))
        struct_err = this;
        return;
      end
      struct_err = struct (...
        "identifier", this.identifier, ...
        "message", this.message, ...
        "stack", this.stack);
    end

  end

  methods (Static)

   function ME = last (reset)

      ## -*- texinfo -*-
      ## @deftypefn {} {@var{ME} =} MException.last ()
      ## @deftypefnx {} {} MException.last (@qcode{"reset"})
      ## Query or reset the last error message information.
      ##
      ## When called without arguments, return a MException object containing
      ## the last error message and other information related to this error.
      ##
      ## If called with the argument "reset", all fields are set to their
      ## default values.
      ## @end deftypefn

      if (nargin)
        if (ischar (reset) && strcmpi (reset, "reset"))
          if (nargout > 0)
            error ("MException: Too many output arguments");
          end
          lasterror ("reset"); # need to reset cause and Correction
        else
          error ("MException: RESET is the only allowed input");
        end
      else
        ME = MException (lasterror ()); # use MException (struct_err) constructor
      end
    end

    function throwWrongNumberInputs()
            ME = MException('FINITELEMENTS:WRONGNUMBERINPUTS', 'The number of input arguments is wrong');
            ME.throwAsCaller();
    end
    function throwWrongNumberOutputs()
            ME = MException('FINITELEMENTS:WRONGNUMBEROUTPUTS', 'The number of output arguments is wrong');
            ME.throwAsCaller();
    end
    function throwWrongFormat()
            ME = MException('FINITELEMENTS:WRONGFORMAT', 'Reading geometry data failed, maybe incorrect format');
            ME.throwAsCaller();
    end
    function throwWrongClass()
            ME = MException('FINITELEMENTS:WRONGCLASS', 'Wrong argument class');
            ME.throwAsCaller();
    end
    function throwWrongSize()
            ME = MException('FINITELEMENTS:WRONGSIZE', 'Wrong sized arguments');
            ME.throwAsCaller();
    end
    function throwEmptyCall()
            ME = MException("MException:EmptyCall",...
                "Constructor must be called at least with one argument.");
            ME.throwAsCaller();
    end
    function throwMustBeVector()
            ME = MException("MException:MustBeVector",...
                "First argument must be a vector.");
            ME.throwAsCaller();
    end
    function throwMustBeMoreThanOnePoint()
            ME = MException("MException:MustBeMoreThanOnePoint",...
                "In [a, b], a must be different from b.");
            ME.throwAsCaller();
    end

    function throwMustBeScalar()
            ME = MException("MException:MustBeScalar",...
                "Second argument must be a positive scalar.");
            ME.throwAsCaller();
    end

    function throwMustBeGreaterThanZero()
            ME = MException("MException:MustBeGreaterThanZero",...
                "Argument hmax must be greater than zero.");
            ME.throwAsCaller();
    end
    function throwTooManyArguments()
            ME = MException("MException:TooManyArguments",...
                "Too many arguments.");
            ME.throwAsCaller();
    end
    % New method for custom error in finiteelements2D
    function throwBilinear3DWrongArgs()
            ME = MException("MException:Bilinear3DWrongArgs",...
                 "The first argument must be a grid2DR object.");
            ME.throwAsCaller();
    end
    function throwWrongArgs()
            ME = MException('finiteElements2D:WRONGARGS', 'The first argument must be a grid2DR object.');
            ME.throwAsCaller();
    end
    % new method for custom error in finiteelements3D
    function throwInvalidCvalSize(cval)
            ME = MException("MException:InvalidCvalSize", ...
                 ["The object cval has " num2str(size(cval, 1)) " Rows."...
                 " If Diffusion Coefficient is scalar, cval must have three equal rows. "...
                 "If it is a diagonal matrix, cval must have these three diagonal elements as rows. "...
                 "If it is a full matrix, cval must have nine rows, which are taken from the matrix columnwise. "...
                 "So if the rows are labeled with c1, c2, c3, c4, and so on, then the matrix entries are"...
                 " [c1 c4 c7; c2 c5 c8; c3 c6 c9]."...
                 " Other syntaxes are not supported. For further help, see aCoefficients() in your used grid class."]);
            ME.throwAsCaller();
    end
    %new method for custom error in grid1D
    function throwWrongFormatGRID1D()
            ME = MException('GRID1D:WRONGFORMAT', 'The format geometry is wrong. Allowed is only double.');
            ME.throwAsCaller();
    end
    function throwRefineMeshError()
            ME = MException('GRID1D:REFINEMESH', 'The number of elements and the index given by the argument do not fit.');
            ME.throwAsCaller();
    end
    function throwWrongNumberOfPoints()
            ME = MException('GRID1D:WRONGNUMBEROFPOINTS', 'The number of grid and data points do not fit.');
            ME.throwAsCaller();
    end
    function throwWrongConvectionFormat()
            ME = MException('GRID1D:WRONGFORMAT', 'The format of convection coefficient is wrong. It must be scalar or a vector of length nElements or nPoints of class double.');
            ME.throwAsCaller();
    end
    function throwWrongConvectionHandleFormat()
            ME = MException('GRID1D:WRONGFORMAT', ...
                'The format of convection coefficient is wrong. Allowed is function_handle, char string or double.');
            ME.throwAsCaller();
    end
    function throwWrongDiffusionFormat()
            ME = MException('GRID1D:WRONGFORMAT', ...
                'The format of diffusion coefficient is wrong. Allowed is function_handle, char string or double.');
            ME.throwAsCaller();
    end
    function throwWrongMassFormatScaler()
            ME = MException('GRID1D:WRONGFORMAT', ...
                'The format of mass coefficient is wrong. It must be scalar or a vector of length nElements or nPoints of class double.');
            ME.throwAsCaller();
    end
    function throwWrongMassFormatHandle()
            ME = MException('GRID1D:WRONGFORMAT', ...
                'The format of mass coefficient is wrong. Allowed is function_handle, char string or double.');
            ME.throwAsCaller();
    end
    function throwWrongSourceFormatScaler()
            ME = MException('GRID1D:WRONGFORMAT', ...
                'The format of source function is wrong. It must be scalar or a vector of length nElements or nPoints of class double.');
            ME.throwAsCaller();
    end
    function throwWrongSourceFormatHandle()
            ME = MException('GRID1D:WRONGFORMAT', ...
                'The format of source function is wrong. Allowed is function_handle, char string, or double.');
            ME.throwAsCaller();
    end
    %new custom methods of grid2D
    function throwWrongEllipseInput()
            ME = MException('ELLIPSE:wrongDefinedInput', ...
                'Error: A and B must be greater than zero.');
            ME.throwAsCaller();
    end
    function throwInvalidGridRefinement()
            ME = MException('GRID2D:REFINEMESH:INVALIDGRID', 'Error when refining the mesh. Check the structure of your mesh. Message was: ');
            ME.throwAsCaller();
    end
    function throwWrongCoefficientDefinitionScaler()
            ME = MException('CCOEFFICIENTS:WRONGCOEFFICIENTDEFINITION', 'b must be a vector');
            ME.throwAsCaller();
    end
    function throwWrongCoefficientDefinition()
            ME = MException('CCOEFFICIENTS:WRONGCOEFFICIENTDEFINITION', 'Wrong coefficient definition');
            ME.throwAsCaller();
    end
    function throwWrongSizeB()
            ME = MException('CCOEFFICIENTS:WRONGSIZE', 'Wrong sized b(1)');
            ME.throwAsCaller();
    end
    function throwWrongFormatB()
            ME = MException('CCOEFFICIENTS:WRONGFORMAT', 'Wrong formated b(1)');
            ME.throwAsCaller();
    end
    function throwNoVector()
            ME = MException('GRID2D:NOVECTOR', 'The argument must be the coordinates of a single point');
            ME.throwAsCaller();
    end
    %new custom methods of grid2DR
    function throwEmptyMesh()
            ME = MException('GRID2DR:EMPTYMESH', 'It makes no sense to refine empty meshes.');
            ME.throwAsCaller();
    end
    function throwWrongElement()
            ME = MException('GRID2D:PLOT:WRONGELEMENT', 'Elements can only be P0, P1, or P2.');
            ME.throwAsCaller();
    end
    function throwCellArrayNotImplemented()
            ME = MException('MATLAB:CellArrayNotImplemented', 'Sorry, Cell array input not yet implemented!');
            ME.throwAsCaller();
    end
    function throwWrongUse3DPR()
            ME = MException('GRID3DPR:WRONGUSE',...
                ['In the 3-argument call, the arguments must ',...
                'be vectors defining the coordinates of the points of the bar.']);
            ME.throwAsCaller();
    end
    %new custom methods of grid3D
    function throwEmptySetCutawayPlot()
            ME = MException('grid3Dpr:cutawayPlot:EmptySet',...
                'The cut is not well defined, check the arguments.');
            ME.throwAsCaller();
    end
    %new custom methods of grid3Dpr
    function throwWrongUse3DPR()
        ME = MException('GRID3DPR:WRONGUSE',...
            ['In the 3 argument call, the arguments must',...
            ' be vectors defining',...
            ' the coordinates of the points of the bar.']);
        ME.throwAsCaller();
    end
    function throwWrongNumberOfArguments3DPR()
        ME = MException('GRID3DPR:WRONGNUMBEROFARGUMENTS',...
            [' The number of arguments must be 3, 6, or 7.',...
            ' See help grid3Dpr.bar.']);
        ME.throwAsCaller();
    end
    function throwACoefficientsNotImplemented()
        ME = MException('grid3Dpr:aCoefficients',...
            ['Not yet implemented, try to use a function',...
            ' for the midpoints of elements']);
        ME.throwAsCaller();
    end
    function throwCellArrayInputNotImplemented()
        ME = MException('MException:CellArrayInputNotImplemented',...
            'Sorry, Cell array input makes here no sense.');
        ME.throwAsCaller();
    end
    % Assuming 'z' is the variable causing the error
    function throwWrongClassError(z)
        ME = MException('MException:WrongClassError',...
            ['The second argument was a ', class(z), '. Allowed is grid1D or double.']);
        ME.throwAsCaller();
    end
    %new custom methods of gridd
    function throwNoGridError()
        ME = MException('MException:NoGridError', 'No grid defined.');
        ME.throwAsCaller();
    end
    function throwNoBoundaryCondition()
            ME = MException('MException:NoBoundaryCondition',...
                'No boundary conditions defined.');
            ME.throwAsCaller();
    end
    function throwWrongBoundaryCondition()
            ME = MException('MException:WrongBoundaryCondition',...
                'The number of boundary conditions in matrix b and in geometry don''t match.');
            ME.throwAsCaller();
    end
    function throwWrongClass()
            ME = MException('MException:WRONGCLASS',...
                'Wrong argument class, check it.');
            ME.throwAsCaller();
    end
    function throwOpNoSense()
            ME = MException('GRIDD:OPNOSENSE',...
                'Size of y is already nElements. Check your code!');
            ME.throwAsCaller();
    end
    function throwOpNoSenseNotN()
            ME = MException('GRIDD:OPNOSENSE',...
                'Size of y is not nPoints.');
            ME.throwAsCaller();
    end
    function throwUnexpectedNumber()
            ME = MException('PDE:point2Center:UnexpectedNumber',...
                 'Unexpected number of points in element');
            ME.throwAsCaller();
    end
    function throwInternalError()
            ME = MException('GRIDD:INTERNAL', 'You found a bug.');
            ME.throwAsCaller();
    end
    function throwOpNoSenseError()
            ME = MException('GRIDD:OPNOSENSE', 'Size of y is not nElements.');
            ME.throwAsCaller();
    end
    function throwIndexErrortoN()
            ME = MException('GRIDD:INDEX', 'Index must be in the range of 1 to nElements.');
            ME.throwAsCaller();
    end
    function throwSenselessActionError()
            ME = MException('GRIDD:SENSELESSACTION', 'You try to rotate an empty grid.');
            ME.throwAsCaller();
    end
    function throwSenselessOperation()
            ME = MException('GRIDD:SENSELESSOPERATION', 'It is not possible to rotate 1D objects.').throwAsCaller;
            ME.throwAsCaller();
    end
    function throwUnknownError()
            ME = MException('GRIDD:UNKNOWNERROR', 'This should never happen.').throwAsCaller;
            ME.throwAsCaller();
    end
    function throwWrongParameter()
            ME = MException('GRIDD:WRONGPARAMETER', 'Scale factor too small or too large.');
            ME.throwAsCaller();
    end
    function throwWrongFormat()
            ME = MException('GRIDD:WRONGFORMAT', 'Input has the wrong format, check it.');
            ME.throwAsCaller();
    end
    function throwNotSupportedDim()
            ME = MException('PDE:NOTSUPPORTEDDIM', 'PDE supports only FE in 1D-3D domains');
            ME.throwAsCaller();
    end
    function throwWrongDataFormat()
            ME = MException('GRIDD:WRONGDATAFORMAT', 'Data y must have nPoints rows.');
            ME.throwAsCaller();
    end
    function throwNotValidFor1D()
            ME = MException('GRIDD:NOTVALID', 'This method makes no sense for 1D grids.');
            ME.throwAsCaller();
    end
    function throwWrongDefinition()
            ME = MException('obj:WRONGDEFINITION', 'You can define only Dirichlet OR Robin BCs.');
            ME.throwAsCaller();
    end
    function throwWrongBoundaryCondition()
            ME = MException('GRIDD:WRONGBOUNDARYCOND', 'The number of boundary conditions in matrix b and in geometry don''t match.');
            ME.throwAsCaller();
    end
    %new custom methods for pde
    function throwForbiddenProperty()
            ME = MException('PDE:FORBIDDEN', 'You cannot set fem property');
            ME.throwAsCaller();
    end
    function throwForbiddenGridProperty()
            ME = MException('PDE:FORBIDDEN', 'You cannot set grid property');
            ME.throwAsCaller();
    end
    function throwMissformedArgument()
            ME = MException('pde:setBoundaryCondition:missformedArgument', 'The argument must be scalar');
            ME.throwAsCaller();
    end
    function throwInvalidArgumentDoubleorCell()
            ME = MException('pde:setBoundaryCondition:missformedArgument',...
                            'The argument must be a double or cell vector of length two.');
            ME.throwAsCaller();
    end
    function throwWrongGradientMode()
            ME = MException('pde:gradient:wrongmode',...
                            'The mode must be *points* or *elements*');
            ME.throwAsCaller();
    end
    function throwWrongGradientFormat()
            ME = MException('pde:gradient:wrongformat',...
                            'Input must be a vector of length *no of points*.');
            ME.throwAsCaller();
    end
    function throwInvalidRefinement()
            ME = MException('PDE:NOTVALID', ['No data in obj.y. pde.refineMesh is for ',...
                            ' refinement of mesh + data.\n',...
                            'Use pde.grid.refineMesh instead']);
            ME.throwAsCaller();
    end
    function throwRefineNotImplemented()
            ME = MException('PDE:NOTIMPLEMENTED', 'Local refinement is still not implemented');
            ME.throwAsCaller();
    end
    function throwWrongNumberArgs()
            ME = MException('PDE:WRONGNUMBERINPUTS', 'Wrong number of arguments.');
            ME.throwAsCaller();
    end
    function throwUnsupportedDimension()
            errorMessage = 'PDE supports only FE in 1D-3D domains.';
            ME = MException('PDE:NOTSUPPORTEDDIM', errorMessage);
            ME.throwAsCaller();
    end
    function throwObjectNotInitialized()
            ME = MException('PDE:NOTINITIALIZED', 'The object is not initialized. Run pde.initialize before calling pde.solve.');
            ME.throwAsCaller();
    end
    function throwWrongDDimension()
            ME = MException('PDE:WRONGD', 'Dimension of D does not fit the dimension of K, M, y, etc.');
            ME.throwAsCaller();
    end
    function throwSolverFailed()
            ME = MException('PDE:SOLVE:SOLVERFAILED', 'Iterative solver failed.');
            ME.throwAsCaller();
    end
    function throwAMGSolverNotFound()
            ME = MException('PDE:AMG:SOLVER_NOT_FOUND', ['Cannot find AMG solver. Be sure that ilupack',...
                ' is installed and configured.']);
            ME.throwAsCaller();
    end
    %new custom methods of plotUtilsTimeDependent
    function throwOnly1DDomain()
            ME = MException('plotUtilsTimeDependent:ONLY1D', 'This method is only defined for 1D domains.');
            ME.throwAsCaller();
    end
    function throwDimensionMismatch()
            ME = MException('PLOTUTILSTIMEDEPENDENT:DIMENSION', 'Time and/or space discretization don''t match data.');
            ME.throwAsCaller();
    end
    function throwWrongThreeDimension()
            ME = MException('PLOTUTILSTIMEDEPENDENT:WRONGDIMENSION', 'Domain must have dimension 1, 2, or 3.');
            ME.throwAsCaller();
    end
    %new custom methods of DoubleT.m geometry 2D
    function throwPositiveScaleFactors()
            ME = MException('GRIDD:POSITIVE', 'Scale factors must be positive.');
            ME.throwAsCaller();
    end
    function throwPositiveHmax()
            ME = MException('GRIDD:POSITIVE', 'hmax must be positive.');
            ME.throwAsCaller();
    end
    %new custom method of Circle.m geometry 2D
    function throwWrongHmax(R, hmax)
            ME = MException('CIRCLE:WRONGHMAX', ['The ratio R/hmax = ', ...
                 num2str(R/hmax), ...
                 ' < 2. Circle is not able to ', ...
                 ' create deformed circles.']);
            ME.throwAsCaller();
    end
    %new custom methods of Bilinear3D of elements
    function throwWrongNumberArgsBilinear3D()
            ME = MException('BILINEAR3D:WRONGNUMBERARGS', 'Wrong number of arguments.');
            ME.throwAsCaller();
    end
    function throwWrongfirstArgs()
            ME = MException('BILINEAR3D:WRONGARGS', 'The first argument must be a grid3Dpr object.');
            ME.throwAsCaller();
    end
    function throwEmptyBoundaryBilinear3D()
            ME = MException('BILINEAR3D:EMPTYBOUNDARY', 'The boundary field is empty, initialize it!');
            ME.throwAsCaller();
    end
    %new custom methods of Lagrange13D of elements
    function throwWrongNumberArgsLagrange13D()
            ME = MException('LAGRANGE13D:WRONGNUMBERARGS', 'Wrong number of arguments.');
            ME.throwAsCaller();
    end
    function throwWrongfirstArgsLagrange13D()
            ME = MException('LAGRANGE13D:WRONGARGS', 'The first argument must be a grid3D object.');
            ME.throwAsCaller();
    end
    %new custom methods for Breakthrough of pdes
    function throwWrongSolverBreakthrough()
            ME = MException('PDE:WRONGSOLVER', 'Non suitable solver for this problem');
            ME.throwAsCaller();
    end
    %new custom methods of Elliptic of pdes
    function throwWrongArgumentClassElliptic()
            ME = MException('ELLIPTIC:WRONGARGUMENTCLASS', 'The argument must be of class char');
            ME.throwAsCaller();
    end
    function throwNotSupportedSolver(solver)
            ME = MException('ELLIPTIC:NOTSUPPORTED', ['The solver ', solver, ' is not supported']);
            ME.throwAsCaller();
    end
    function throwNoCoefficientFunctionsSet()
            ME = MException('PDE:COEFFIENCETS', 'No coefficient functions set.');
            ME.throwAsCaller();
    end
    function throwWrongNumberArgumentsEllipticAdaptive()
            ME = MException('PDE:WRONGNUMBERARGUMENTS', 'The number of arguments must be zero or three.');
            ME.throwAsCaller();
    end
    function throwNotInitializedEllipticAdaptive()
            ME = MException('PDE:NOTINITIALIZED', 'The object is not initialized.');
            ME.throwAsCaller();
    end
    function throwWrongClassForAdaption()
            ME = MException('PDE:WRONGCLASS', 'For adaption, all data must be function_handles or scalar.');
            ME.throwAsCaller();
    end
    function throwWrongClassForAdaptiondouble()
            ME = MException('PDE:WRONGCLASS', 'For adaption, all data must be function_handles or double.');
            ME.throwAsCaller();
    end
    %new custom method of embededGridObject of postprocessing
    function throwWrongNumberArgumentsEmbeddedGridObject()
            ME = MException('EMBEDEDGRIDOBJECT:WRONGNUMBERARGUMENTS', 'This call needs at least two arguments.');
            ME.throwAsCaller();
    end
  end
end

%!test
%! ME = MException ("octave:error", "error message");
%! ME = MException ("octave:function:error", "error message");
%! ME = MException ("octave:error", "file '%s'", "data.dat");
%! identifier = ME.identifier;
%! message = ME.message;
%! stack = ME.stack;
%! cause = ME.cause;
%! Correction = ME.Correction;

## Test error throwing
%!error
%! ME = MException ("octave:error", "error message");
%! throw (ME);
%!error
%! ME = MException ("octave:error", "error message");
%! throwAsCaller (ME);
%!error
%! try
%!   error ("error message");
%! catch ME
%!   rethrow (ME);
%! end_try_catch

## Test MException.addCause ()
%!test
%! ME = MException ("octave:error", "error message");
%! ME = ME.addCause (ME);

## Test MException.addCorrection ()
%!test
%! ME = MException ("octave:error", "error message");
%! %C  = matlab.lang.correction.Correction ();
%! %ME = ME.addCorrection (C);

## Test MException.getReport ()
%!test
%! ME = MException ("octave:error", "error message");
%! report = ME.getReport ();
%! assert (ischar (report));
%! report = ME.getReport ("basic");
%! assert (report, ME.message);
%! report = ME.getReport ("extended");
%! assert (ischar (report));
%! report = ME.getReport ("extended", "hyperlinks", "on");
%! report = ME.getReport ("extended", "hyperlinks", "off");
%! report = ME.getReport ("extended", "hyperlinks", "default");

## Test MException.last ()
%!test
%! ME = MException.last ();
%! M = ME.last ();
%! MException.last ("reset");

## Test input validation
%!error MException;
%!error MException ("octave:error");
%!error MException ("octave", "error message");
%!error MException ("octave:err!", "error message");
%!error MException ("octave:1st", "error message");
%!error MException ("1st:error", "error message");
%!error MException (1, "error message");
%!error MException ("octave:error", NaN);
%!error MException ("octave:error", "error message").getReport ("full");
%!error MException ("octave:error", "error message").getReport ("basic", "hyperlinks");
%!error MException ("octave:error", "error message").getReport ("basic", "a", "b");
%!error MException.last ("random");
%!error ME = MException.last ("reset");
%!error
%! ME = MException ("octave:error", "error message");
%! ME.identifier = "octave:error";
%!error
%! ME = MException ("octave:error", "error message");
%! ME.message = "error message";
%!error
%! ME = MException ("octave:error", "error message");
%! ME.stack = struct ();
%!error
%! ME = MException ("octave:error", "error message");
%! ME.cause = [];
%!error
%! ME = MException ("octave:error", "error message");
%! ME.Correction = [];

