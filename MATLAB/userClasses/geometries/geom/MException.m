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
    wrongNumberInputs;
    wrongNumberOutputs;
    wrongFormat;
    wrongClass;
    wrongSize;
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

    function instance = getWrongNumberInputs()
      instance = MException('FINITELEMENTS:WRONGNUMBERINPUTS', 'The number of input arguments is wrong');
    end

    function instance = getWrongNumberOutputs()
      instance = MException('FINITELEMENTS:WRONGNUMBEROUTPUTS', 'The number of output arguments is wrong');
    end

    function instance = getWrongFormat()
      instance = MException('FINITELEMENTS:WRONGFORMAT', 'Reading geometry data failed, maybe incorrect format');
    end

    function instance = getWrongClass()
      instance = MException('FINITELEMENTS:WRONGCLASS', 'Wrong argument class');
    end

    function instance = getWrongSize()
      instance = MException('FINITELEMENTS:WRONGSIZE', 'Wrong sized arguments');
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

