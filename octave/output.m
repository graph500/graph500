function output (SCALE, edgefactor, NBFS, kernel_1_time, kernel_2_time, kernel_2_nedge)
  printf ("SCALE: %d\n", SCALE);
  printf ("edgefactor: %d\n", edgefactor);
  printf ("NBFS: %d\n", NBFS);
  printf ("construction_time: %20.17e\n", kernel_1_time);

  S = statistics (kernel_2_time);
  printf ("min_time: %20.17e\n", S(1));
  printf ("firstquartile_time: %20.17e\n", S(2));
  printf ("median_time: %20.17e\n", S(3));
  printf ("thirdquartile_time: %20.17e\n", S(4));
  printf ("max_time: %20.17e\n", S(5));
  printf ("mean_time: %20.17e\n", S(6));
  printf ("stddev_time: %20.17e\n", S(7));

  S = statistics (kernel_2_nedge);
  printf ("min_nedge: %20.17e\n", S(1));
  printf ("firstquartile_nedge: %20.17e\n", S(2));
  printf ("median_nedge: %20.17e\n", S(3));
  printf ("thirdquartile_nedge: %20.17e\n", S(4));
  printf ("max_nedge: %20.17e\n", S(5));
  printf ("mean_nedge: %20.17e\n", S(6));
  printf ("stddev_nedge: %20.17e\n", S(7));

  TEPS = kernel_2_nedge ./ kernel_2_time;
  N = length (TEPS);
  S = statistics (TEPS);
  S(6) = mean (TEPS, 'h');
  %% Harmonic standard deviation from:
  %% Nilan Norris, The Standard Errors of the Geometric and Harmonic
  %% Means and Their Application to Index Numbers, 1940.
  %% http://www.jstor.org/stable/2235723
  tmp = zeros (N, 1);
  tmp(TEPS > 0) = 1./TEPS(TEPS > 0);
  tmp = tmp - 1/S(6);
  S(7) = (sqrt (sum (tmp.^2)) / (N-1)) * S(6)^2;
  
  printf ("min_TEPS: %20.17e\n", S(1));
  printf ("firstquartile_TEPS: %20.17e\n", S(2));
  printf ("median_TEPS: %20.17e\n", S(3));
  printf ("thirdquartile_TEPS: %20.17e\n", S(4));
  printf ("max_TEPS: %20.17e\n", S(5));
  printf ("harmonic_mean_TEPS: %20.17e\n", S(6));
  printf ("harmonic_stddev_TEPS: %20.17e\n", S(7));
