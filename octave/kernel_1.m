function G = kernel_1 (ij)
%% Compute a sparse adjacency matrix representation
%% of the graph with edges from ij.

  %% Remove self-edges.
  ij(:, ij(1,:) == ij(2,:)) = [];
  %% Adjust away from zero labels.
  ij = ij + 1;
  %% Find the maximum label for sizing.
  N = max (max (ij));
  %% Create the matrix, ensuring it is square.
  G = sparse (ij(1,:), ij(2,:), ones (1, size (ij, 2)), N, N);
  %% Symmetrize to model an undirected graph.
  G = spones (G + G.');
