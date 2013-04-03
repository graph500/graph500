function G = kernel_1 (ijw)
## Compute a sparse adjacency matrix representation
## of the graph with edges from ij.

  ## Remove self-edges.
  ijw(ijw(:, 1) == ijw(:, 2), :) = [];
  ## Adjust away from zero labels.
  ijw(:, [1 2]) = ijw(:, [1 2]) + 1;
  ## Find the maximum label for sizing.
  N = max (max (ijw(:, [1 2])));
  ## Order into a single triangle.
  mask = ijw(:, 1) < ijw(:, 2);
  ijw(mask, [1 2]) = ijw(mask, [2 1]);
  ## Create the matrix, ensuring it is square.
  G = sparse (ijw(:, 1), ijw(:, 2), ijw(:, 3), N, N);
  ## Symmetrize to model an undirected graph.
  G = G + G.';
endfunction
