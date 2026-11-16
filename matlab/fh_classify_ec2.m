function cls = fh_classify_ec2(fR1k, fR3k)
%FH_CLASSIFY_EC2 Class of SFRC per EN 1992-1-1:2023 Table L.2.
series = [1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 7.0 8.0];
N = max(series(series <= fR1k));
if isempty(N), cls = 'n.c.'; return; end
lim = [0.5 0.7 0.9 1.1 1.3]; let = 'abcde';
k = find(fR3k >= lim*N, 1, 'last');
if isempty(k), cls = sprintf('%g -', N); else, cls = sprintf('%g %s', N, let(k)); end
end
