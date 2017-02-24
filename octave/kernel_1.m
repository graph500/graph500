function G = kernel_1 (ijw)
%% Compute a sparse adjacency matrix representation
%% of the graph with edges from ijw.

  %% Remove self-edges.
  ijw(:, ijw(1,:) == ijw(2,:)) = [];
  %% Adjust away from zero labels.
  ijw(1:2,:) = ijw(1:2,:) + 1;
  %% Order into a single triangle
  mask = ijw(1, :) < ijw(2, :);
  ijw([1 2], mask) = ijw([2 1], mask);
  %% Find the maximum label for sizing.
  N = max (max (ijw(1:2,:)));
  %% Create the matrix, ensuring it is square.
  G = accumarray ([ijw(1,:)', ijw(2,:)'], ijw(3,:)', [N, N], @min, 0, true);
  %% Symmetrize to model an undirected graph.
  G = G + G.';

