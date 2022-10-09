function S = RenameField(S, Old, New)
% RenameField - Rename a field of a struct
% T = RenameField(S, Old, New)
% INPUT:
%   S: Struct or struct array.
%   Old, New: CHAR vectors or cell string. In the M-Version of this function
%      strings are allowed also. The names must be valid Matlab symbols:
%      <= 63 characters, first character is a letter, others characters are
%      alpha-numeric or the underscore.
%      If a name in Old exist in S, it is renamed to the corresponding name in
%      New. Not existing names are ignored.
% OUTPUT:
%   T: Struct S with renamed fields.
%
% EXAMPLES:
%   S.A = 1; S.B = 2;
%   T = RenameField(S, 'B', 'C');  % T.A = 1, T.C = 2
%
% NOTES:
% * The C-Mex-function is much faster, but does not accept strings, only CHAR
%   vectors and cell string.
% * Hardcore programmers can omit the validity checks of the new name in the
%   C-Mex function: All names up to 63 ASCII characters are allowed, even '*',
%   ' ' and the empty char vector ''. This does not crash Matlab in my
%   experiments. In R2018b such fields can be accessed with dynamic fieldnames:
%   S.(''), S.('*') works!
% * The check for multiple field names can be omitted also. Dynamic field names
%   pick the first occurence in this case, but e.g. RMFIELD stops with an error.
%   Other functions might crash.
% * This function was created after a discussion in Loren's blog:
%   http://blogs.mathworks.com/
%          loren/2010/05/13/rename-a-field-in-a-structure-array
%
% COMPILATION:
%   InstallMex('RenameField', 'uTest_RenameField')
% See RenameField.c for more details.
%
% Tested: Matlab 2009a, 2015b(32/64), 2016b, 2018b, Win7/10
% Author: Jan Simon, Heidelberg, (C) 2006-2022 matlab.2010(a)n(MINUS)simon.de
%
% See also CELL2STRUCT, STRUCT, GENVARNAME, RMFIELD.

% $JRev: R-l V:011 Sum:A/2XtAQ6eCju Date:09-Jun-2022 01:08:15 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLStruct\RenameField.m $
% History:
% 001: 19-Aug-2010 00:11, Created after discussion in Loren's blog.
% 002: 08-Jun-2022 22:52, STRING type included in M-verion.

% Initialize: ==================================================================
% Do the work: =================================================================

% This is an implementation as M-code. Prefer the mex file, which is 50% (S has
% 1 field only) to 95% (S has 1000 fields) faster.

% Comment this out, if you want to use the M-version:
% error(['JSimon:', mfilename, ':NoMex'], ...
%       ['*** ', mfilename, ': Cannot find compiled Mex file!']);

% Under Matlab 6.5 CELL2STRUCT accepted names with more than 63 characters. But
% the later recognition fails!

if isempty(S) && isa(S, 'double')  % Accept [] as empty struct without fields
   return;
end

Data  = struct2cell(S);
Field = fieldnames(S);
if ischar(Old)
   Field(strcmp(Field, Old)) = {New};

elseif iscellstr(Old)   %#ok<ISCLSTR>
   for iField = 1:numel(Old)
      match = strcmp(Field, Old{iField});
      if any(match)
         Field{match} = New{iField};
      end
   end
   
elseif isa(Old, 'string')
   for iField = 1:numel(Old)
      match = strcmp(Field, Old(iField));
      if any(match)
         Field{match} = char(New(iField));
      end
   end
   
else
   error(['JSimon:', mfilename, ':BadInputType'], ...
      '*** %s: Names must be CHAR vectors, cell strings or strings!', ...
      mfilename);
end

S = cell2struct(Data, Field);

end