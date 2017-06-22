function output (SCALE, NBFS, NSSSP, kernel_1_time, kernel_2_time, kernel_2_nedge, kernel_3_time, kernel_3_nedge)
  printf ("SCALE: %d\n", SCALE);
  printf ("NBFS: %d\n", NBFS);
  printf ("construction_time: %20.17e\n", kernel_1_time);

  S = statistics (kernel_2_time);
  printf ("bfs_min_time: %20.17e\n", S(1));
  printf ("bfs_firstquartile_time: %20.17e\n", S(2));
  printf ("bfs_median_time: %20.17e\n", S(3));
  printf ("bfs_thirdquartile_time: %20.17e\n", S(4));
  printf ("bfs_max_time: %20.17e\n", S(5));
  printf ("bfs_mean_time: %20.17e\n", S(6));
  printf ("bfs_stddev_time: %20.17e\n", S(7));

  S = statistics (kernel_2_nedge);
  printf ("bfs_min_nedge: %20.17e\n", S(1));
  printf ("bfs_firstquartile_nedge: %20.17e\n", S(2));
  printf ("bfs_median_nedge: %20.17e\n", S(3));
  printf ("bfs_thirdquartile_nedge: %20.17e\n", S(4));
  printf ("bfs_max_nedge: %20.17e\n", S(5));
  printf ("bfs_mean_nedge: %20.17e\n", S(6));
  printf ("bfs_stddev_nedge: %20.17e\n", S(7));

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
  
  printf ("bfs_min_TEPS: %20.17e\n", S(1));
  printf ("bfs_firstquartile_TEPS: %20.17e\n", S(2));
  printf ("bfs_median_TEPS: %20.17e\n", S(3));
  printf ("bfs_thirdquartile_TEPS: %20.17e\n", S(4));
  printf ("bfs_max_TEPS: %20.17e\n", S(5));
  printf ("bfs_harmonic_mean_TEPS: %20.17e\n", S(6));
  printf ("bfs_harmonic_stddev_TEPS: %20.17e\n", S(7));

  S = statistics (kernel_3_time);
  printf ("sssp_min_time: %20.17e\n", S(1));
  printf ("sssp_firstquartile_time: %20.17e\n", S(2));
  printf ("sssp_median_time: %20.17e\n", S(3));
  printf ("sssp_thirdquartile_time: %20.17e\n", S(4));
  printf ("sssp_max_time: %20.17e\n", S(5));
  printf ("sssp_mean_time: %20.17e\n", S(6));
  printf ("sssp_stddev_time: %20.17e\n", S(7));

  S = statistics (kernel_3_nedge);
  printf ("sssp_min_nedge: %20.17e\n", S(1));
  printf ("sssp_firstquartile_nedge: %20.17e\n", S(2));
  printf ("sssp_median_nedge: %20.17e\n", S(3));
  printf ("sssp_thirdquartile_nedge: %20.17e\n", S(4));
  printf ("sssp_max_nedge: %20.17e\n", S(5));
  printf ("sssp_mean_nedge: %20.17e\n", S(6));
  printf ("sssp_stddev_nedge: %20.17e\n", S(7));

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
  
  printf ("sssp_min_TEPS: %20.17e\n", S(1));
  printf ("sssp_firstquartile_TEPS: %20.17e\n", S(2));
  printf ("sssp_median_TEPS: %20.17e\n", S(3));
  printf ("sssp_thirdquartile_TEPS: %20.17e\n", S(4));
  printf ("sssp_max_TEPS: %20.17e\n", S(5));
  printf ("sssp_harmonic_mean_TEPS: %20.17e\n", S(6));
  printf ("sssp_harmonic_stddev_TEPS: %20.17e\n", S(7));
