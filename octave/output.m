function output (machine, SCALE, edgefactor,
		 root, kernel_1_time,
		 kernel_2_time, kernel_2_dmax,
		 kernel_2_verify_time,
		 kernel_3_time, kernel_3_dmax,
		 kernel_3_verify_time,
		 comment=[])
  printf ("MACHINE: %s\n", machine);
  if !isempty (comment), printf ("COMMENT: %s\n", comment); endif
  printf ("IMPLEMENTATION: %s\n",
	  "Pseudo-reference, unoptimized Octave");
  printf ("SCALE: %d\n", SCALE);
  printf ("EDGEFACTOR: %d\n", edgefactor);
  printf ("NROOT: %d\n", length (root));
  printf ("PRNGCHECK: %d\n", PRNGCHECK (SCALE, edgefactor));
  printf ("K1TIME: %11.8e\n", kernel_1_time);

  NV = 2**SCALE;
  NE = edgefactor * NV;

  ## Extra, not required fields...
  [mn2, sd2] = avg_teps (NE, kernel_2_time);
  printf ("K2TEPSMEAN: %11.8e\n", mn2);
  printf ("K2TEPSSTDDEV: %11.8e\n", sd2);
  [mn3, sd3] = avg_teps (NE, kernel_3_time);
  printf ("K3TEPSMEAN: %11.8e\n", mn3);
  printf ("K3TEPSSTDDEV: %11.8e\n", sd3);

  printf ("\nroot,k2time,k2max,k2vtime,k3time,k3max,k3vtime\n");
  for k=1:length (root),
    printf ("%d,%11.8e,%d,%11.8e,%11.8e,%d,%11.8e\n",
	    root(k),
	    kernel_2_time(k), kernel_2_dmax(k),
	    kernel_2_verify_time(k),
	    kernel_3_time(k), kernel_3_dmax(k),
	    kernel_3_verify_time(k));
  endfor
  ## The verification times are not required but can be
  ## informative.
endfunction

function [mn, sd] = avg_teps (NE, time)
  TEPS = NE ./ time;
  mn = mean (TEPS, 'h');
  N = length (time);
  ## Harmonic standard deviation from:
  ## Nilan Norris, The Standard Errors of the Geometric and Harmonic
  ## Means and Their Application to Index Numbers, 1940.
  ## http://www.jstor.org/stable/2235723
  tmp = zeros (N, 1);
  tmp(TEPS > 0) = 1./TEPS(TEPS > 0);
  tmp = tmp - 1/mn;
  sd = (sqrt (sum (tmp.^2)) / (N-1)) * mn^2;
endfunction
