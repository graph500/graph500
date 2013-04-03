function [v1, v2, w] = rmat_edge (k, SCALE)
  ## Set initiator probabilities.
  [A, B] = deal (0.55, 0.1);
  ## Noise factor for perturbing the initiator.
  noisefact = 0.1;

  ## Collect all the PRNG outputs used by the
  ## vertices at k.  Each call returns four
  ## pseudo-random numbers used alternately for
  ## the parameter perturbation and the quadrant
  ## over *two* scales.  If SCALE is odd, this
  ## will waste two generated numbers.
  rnd = zeros (1, 2*(SCALE + mod (SCALE, 2)));
  for scl=0:2:(SCALE-1),
    idx = 1 + ((4*floor ((scl+1)/2)):(4*floor ((scl+2)/2)-1));
    rnd(idx) = PRNG (k, 1+floor (scl/2));
  endfor
  rnd = reshape (rnd(1:2*SCALE), 2, SCALE);
  if mod (SCALE, 2), rnd(:,SCALE+1) = []; endif
  rnd = rnd.'; # Silly optimization for
  # column-major ordering.

  mu = noisefact * (2 * rnd(:, 1) - 1);
  As = A * (1 - 2 * mu / (1 - 2*B));
  Bs = B * (1 + mu);

  ## Cast the darts into quadrants using the
  ## perturbed parameters.
  scl = 2.^(0:SCALE-1).';
  v1 = sum ((rnd(:, 2) >= As + Bs) .* scl);
  v2 = sum ((or (and (rnd(:, 2) >= As,
		      rnd(:, 2) < As + Bs),
		 rnd(:, 2) >= As + 2*Bs)) .* scl);
endfunction
