function opts = sidParseOptions(defs, args)
% SIDPARSEOPTIONS Parse name-value pairs from varargin.
%
%   opts = sidParseOptions(defs, args)
%
%   Parses a cell array of name-value pairs and returns a struct with
%   the parsed values. Unknown option names produce an error. Option
%   name matching is case-insensitive.
%
%   INPUTS:
%     defs - Struct whose field names define the valid options and whose
%            field values provide the defaults.
%     args - Cell array of name-value pairs (typically varargin from the
%            calling function).
%
%   OUTPUTS:
%     opts - Struct with the same fields as defs, updated with any
%            values supplied in args.
%
%   EXAMPLES:
%     defs.WindowSize = 30;
%     defs.Frequencies = [];
%     defs.SampleTime = 1.0;
%     opts = sidParseOptions(defs, {'WindowSize', 50, 'SampleTime', 0.01});
%     % opts.WindowSize == 50, opts.Frequencies == [], opts.SampleTime == 0.01
%
%   SPECIFICATION:
%     Unified option parsing for cross-language portability. Maps directly
%     to Python **kwargs and Julia keyword arguments.
%
%   See also: sidValidateData
%
%   Changelog:
%   2026-04-04: First version by Pedro Lourenço.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenço, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    opts = defs;
    names = fieldnames(defs);
    namesLower = lower(names);

    k = 1;
    while k <= length(args)
        if ~ischar(args{k})
            error('sid:badInput', ...
                'Expected option name (string) at position %d.', k);
        end
        if k + 1 > length(args)
            error('sid:badInput', ...
                'Option ''%s'' has no corresponding value.', args{k});
        end

        idx = find(strcmp(lower(args{k}), namesLower), 1);
        if isempty(idx)
            error('sid:unknownOption', 'Unknown option: %s', args{k});
        end

        opts.(names{idx}) = args{k + 1};
        k = k + 2;
    end
end
