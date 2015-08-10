function output (SCALE, NBFS, NSSSP, kernel_1_time, kernel_2_time, kernel_2_nedge, kernel_3_time, kernel_3_nedge)
  printf ("SCALE: %d\n", SCALE);
  printf ("NBFS: %d\n", NBFS);
  printf ("NSSSP: %d\n", NSSSP);
  printf ("construction_time: %20.17e\n", kernel_1_time);

  S = statistics (kernel_2_time);
  printf ("min_k2time: %20.17e\n", S(1));
  printf ("firstquartile_k2time: %20.17e\n", S(2));
  printf ("median_k2time: %20.17e\n", S(3));
  printf ("thirdquartile_k2time: %20.17e\n", S(4));
  printf ("max_k2time: %20.17e\n", S(5));
  printf ("mean_k2time: %20.17e\n", S(6));
  printf ("stddev_k2time: %20.17e\n", S(7));

  S = statistics (kernel_2_nedge);
  printf ("min_k2nedge: %20.17e\n", S(1));
  printf ("firstquartile_k2nedge: %20.17e\n", S(2));
  printf ("median_k2nedge: %20.17e\n", S(3));
  printf ("thirdquartile_k2nedge: %20.17e\n", S(4));
  printf ("max_k2nedge: %20.17e\n", S(5));
  printf ("mean_k2nedge: %20.17e\n", S(6));
  printf ("stddev_k2nedge: %20.17e\n", S(7));

  K2TEPS = kernel_2_nedge ./ kernel_2_time;
  K2N = length (K2TEPS);
  S = statistics (K2TEPS);
  S(6) = mean (K2TEPS, 'h');
  %% Harmonic standard deviation from:
  %% Nilan Norris, The Standard Errors of the Geometric and Harmonic
  %% Means and Their Application to Index Numbers, 1940.
  %% http://www.jstor.org/stable/2235723
  k2tmp = zeros (K2N, 1);
  k2tmp(K2TEPS > 0) = 1./K2TEPS(K2TEPS > 0);
  k2tmp = k2tmp - 1/S(6);
  S(7) = (sqrt (sum (k2tmp.^2)) / (K2N-1)) * S(6)^2;
  
  printf ("min_K2TEPS: %20.17e\n", S(1));
  printf ("firstquartile_K2TEPS: %20.17e\n", S(2));
  printf ("median_K2TEPS: %20.17e\n", S(3));
  printf ("thirdquartile_K2TEPS: %20.17e\n", S(4));
  printf ("max_K2TEPS: %20.17e\n", S(5));
  printf ("harmonic_mean_K2TEPS: %20.17e\n", S(6));
  printf ("harmonic_stddev_K2TEPS: %20.17e\n", S(7));

  S = statistics (kernel_3_time);
  printf ("min_k3time: %20.17e\n", S(1));
  printf ("firstquartile_k3time: %20.17e\n", S(2));
  printf ("median_k3time: %20.17e\n", S(3));
  printf ("thirdquartile_k3time: %20.17e\n", S(4));
  printf ("max_k3time: %20.17e\n", S(5));
  printf ("mean_k3time: %20.17e\n", S(6));
  printf ("stddev_k3time: %20.17e\n", S(7));

  S = statistics (kernel_3_nedge);
  printf ("min_k3nedge: %20.17e\n", S(1));
  printf ("firstquartile_k3nedge: %20.17e\n", S(2));
  printf ("median_k3nedge: %20.17e\n", S(3));
  printf ("thirdquartile_k3nedge: %20.17e\n", S(4));
  printf ("max_k3nedge: %20.17e\n", S(5));
  printf ("mean_k3nedge: %20.17e\n", S(6));
  printf ("stddev_k3nedge: %20.17e\n", S(7));

  K3TEPS = kernel_3_nedge ./ kernel_3_time;
  K3N = length (K3TEPS);
  S = statistics (K3TEPS);
  S(6) = mean (K3TEPS, 'h');
  %% Harmonic standard deviation from:
  %% Nilan Norris, The Standard Errors of the Geometric and Harmonic
  %% Means and Their Application to Index Numbers, 1940.
  %% http://www.jstor.org/stable/2235723
  k3tmp = zeros (K3N, 1);
  k3tmp(K3TEPS > 0) = 1./K3TEPS(K3TEPS > 0);
  k3tmp = k3tmp - 1/S(6);
  S(7) = (sqrt (sum (k3tmp.^2)) / (K3N-1)) * S(6)^2;
  
  printf ("min_K3TEPS: %20.17e\n", S(1));
  printf ("firstquartile_K3TEPS: %20.17e\n", S(2));
  printf ("median_K3TEPS: %20.17e\n", S(3));
  printf ("thirdquartile_K3TEPS: %20.17e\n", S(4));
  printf ("max_K3TEPS: %20.17e\n", S(5));
  printf ("harmonic_mean_K3TEPS: %20.17e\n", S(6));
  printf ("harmonic_stddev_K3TEPS: %20.17e\n", S(7));
