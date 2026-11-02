function S = read_yaml(filename)
%READ_YAML  minimal YAML reader (MATLAB / GNU Octave, no toolboxes)
%
%   S = READ_YAML(FILENAME) parses the subset of YAML used by the config
%   files of this folder into a nested struct:
%     - block mappings defined by indentation ("key: value" / "key:")
%     - scalars: numbers (incl. 1.0e-8), true/false, null, quoted or plain
%       strings
%     - flow sequences of scalars: [a, b, c] (numeric -> row vector,
%       otherwise cell array of strings)
%     - comments starting with '#' (outside quotes) and blank lines
%   Anchors, multi-line strings and block sequences ("- item") are not
%   supported and raise an error.

txt = fileread(filename);
lines = regexp(txt, '\r?\n', 'split');

S      = struct();
stack  = {};        % path of keys of the currently open mappings
indent = -1;        % indentation of each open level
for n = 1:numel(lines)
    raw  = strip_comment(lines{n});
    if isempty(strtrim(raw)), continue; end
    ind  = numel(raw) - numel(regexprep(raw, '^ +', ''));
    if any(raw == sprintf('\t'))
        error('read_yaml: tabs are not allowed (%s:%d)', filename, n);
    end
    body = strtrim(raw);
    if body(1) == '-'
        error('read_yaml: block sequences are not supported (%s:%d)', filename, n);
    end
    tok = regexp(body, '^([A-Za-z_][\w\-]*)\s*:(.*)$', 'tokens', 'once');
    if isempty(tok)
        error('read_yaml: cannot parse line %d of %s: "%s"', n, filename, body);
    end
    key = strrep(tok{1}, '-', '_');
    val = strtrim(tok{2});

    while ind <= indent(end)                 % close deeper levels
        indent(end) = [];
        stack(end)  = [];
    end
    if isempty(val)                          % opens a new mapping
        S = setfield(S, stack{:}, key, struct());
        stack{end+1}  = key;  %#ok<AGROW>
        indent(end+1) = ind;  %#ok<AGROW>
    else
        S = setfield(S, stack{:}, key, parse_value(val));
    end
end
end

function s = strip_comment(s)
q = '';
for i = 1:numel(s)
    c = s(i);
    if isempty(q) && (c == '"' || c == '''')
        q = c;
    elseif ~isempty(q) && c == q
        q = '';
    elseif isempty(q) && c == '#' && (i == 1 || any(s(i-1) == ' '))
        s = s(1:i-1);
        return;
    end
end
end

function v = parse_value(s)
s = strtrim(s);
if numel(s) >= 2 && s(1) == '[' && s(end) == ']'
    inner = strtrim(s(2:end-1));
    if isempty(inner), v = []; return; end
    parts = strtrim(strsplit(inner, ','));
    vals  = cellfun(@parse_value, parts, 'UniformOutput', false);
    if all(cellfun(@(x) isnumeric(x) && isscalar(x), vals))
        v = cell2mat(vals);
    else
        v = vals;
    end
    return;
end
if numel(s) >= 2 && (s(1) == '"' || s(1) == '''') && s(end) == s(1)
    v = s(2:end-1);
    return;
end
switch lower(s)
    case {'true', 'yes', 'on'},   v = true;   return;
    case {'false', 'no', 'off'},  v = false;  return;
    case {'null', '~'},           v = [];     return;
end
if ~isempty(regexp(s, '^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$', 'once'))
    v = str2double(s);
else
    v = s;
end
end
