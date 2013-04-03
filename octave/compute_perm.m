function [Z, Zinv] = compute_perm (NE)
  for Z = (1+floor(3*NE/4)):(NE-1),
    [g, Zinv] = gcd (Z, NE);
    if 1 == g,
      assert (1 == mod (Z * Zinv, NE));
      return;
    endif
  endfor
endfunction
