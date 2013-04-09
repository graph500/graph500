function root = sample_roots (NV, NROOT, NE)
  NROOT = min (NROOT, NV);
  NE = int64 (NE);

  root = -ones (1, NROOT);

  ## Method A in Jeffrey Scott Vitter, "An
  ## Efficient Algorithm for Sequential Random
  ## Sampling," ACM Transactions on Mathematical
  ## Software, 13(1), March 1987, 58-67.
  N = NV;
  top = NV - NROOT;
  m = 1;
  cur = 0;
  for m=1:NROOT-1,
    rv = dpPRNG (NE, m-1);
    r = rv(1);
    S = 0;
    quot = top / N;
    while quot > r,
      S += 1;
      top -= 1;
      N -= 1;
      quot *= top / N;
    endwhile
    cur += S+1;
    root(m) = cur;
    N -= 1;
  endfor
  rv = dpPRNG (NE, NROOT-1);
  r = rv(1);
  S = floor (N * r);
  cur += S+1;
  root(NROOT) = cur;
  root -= 1; # Zero-indexed.
  assert (root >= 0 && root < NV);
endfunction
