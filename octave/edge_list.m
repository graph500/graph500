function ijw = edge_list (ne_begin, ne_len,
			  SCALE, NE, maxweight)

  NV = 2**SCALE;
  ijw = zeros (3, ne_len);
  [Z, Zinv] = compute_perm (NE);

  for t = 1:ne_len,
    ## kp = k', location in the global edge list.
    kp = ne_begin + t - 1;
    k = mod (Zinv * kp, NE);

    ## Generate four pseudo-random numbers, but use
    ## only one for the weight.
    rnd = PRNG (k, 0);
    w = ceil (rnd(1) * maxweight);

    if k < NV,
      [v1, v2] = tree_edge (k);
    else
      [v1, v2] = rmat_edge (k, SCALE);
    endif
    v1 = scramble (v1, SCALE);
    v2 = scramble (v2, SCALE);
    ijw(:, t) = [v1; v2; w]; # Location kp "globally."
  endfor
  ijw = ijw.'; # Switch to columns for i, j, w.

endfunction
