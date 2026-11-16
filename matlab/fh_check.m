function ok = fh_check(name, val, ref, tol)
%FH_CHECK Prints a PASS/FAIL line comparing val with ref (relative tolerance tol).
ok = abs(val - ref) <= tol*abs(ref);
if ok, st = 'PASS'; else, st = 'FAIL'; end
fprintf('%-38s %10.4f  ref %10.4f  %s\n', name, val, ref, st);
end
