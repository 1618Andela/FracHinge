function s = fh_read_numbers(fname)
  txt = fileread(fname);
  if exist('jsondecode', 'builtin') || exist('jsondecode', 'file')
    s = jsondecode(txt);
  else
    s = struct(); toks = regexp(txt, '"(\w+)":\s*"([^"]*)"', 'tokens');
    for k = 1:numel(toks), s.(toks{k}{1}) = toks{k}{2}; end
  end
end
