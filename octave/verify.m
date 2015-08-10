function out = verify (SCALE, parent, ijw, root, d, is_bfs)
  out = 1;

  ## Adjust to 1-offset.
  parent = parent + 1;
  root = root + 1;
  ijw(:, [1 2]) = ijw(:, [1 2]) + 1;
  ## Remove self-loops.
  ijw(ijw(:, 1) == ijw(:, 2), :) = [];

  N = max (max (ijw(:, [1 2])));

  if parent(root) != root || \
     sum ((1:N).' == parent) != 1,
    ## There is not a unique root.
    out = 0;
    return;
  end

  if size (parent, 1) != N || \
     any (parent <= 0 || parent > N),
    ## Not every vertex is included, or parent out
    ## of range.
    out = -1;
    return;
  end

  if nargin < 4,
    ## Compute the depth vector.
    d = ones (size (parent));
    d(root) = 0;
    P = parent;
    slice = find (P != root);
    while !isempty (slice),
      d(slice) += 1;
      P(slice) = P(P(slice));
      slice = slice(find (P(slice) != root));
      if any (d > N),
        ## There must be a cycle in the tree.
        out = -2;
        return;
      end
    end
  end
  if nargin < 5,
    ## Assume we're verifying BFS.
    is_bfs = 1;
  endif

  ## Order vertex tuples to point away from the
  ## root.
  mask = d(ijw(:, 1)) > d(ijw(:, 2));
  ijw(mask, [1, 2]) = ijw(mask, [2, 1]);
  assert (d(ijw(:, 1)) <= d(ijw(:, 2)));

  tree_nodes = unique (parent);
  parent_child_edge_list = \
    find(ijw(:, 1) == parent(ijw(:, 2)));

  ## Check that every root-facing vertex in a
  ## parent-child edge is an internal tree node.
  if length (unique
	       (ijw(parent_child_edge_list, 1))) \
    != length (tree_nodes),
    ## Some tree edges are not in the input edge
    ## list.
    out = -3;
  endif

  ## Coping with duplicate edges without collapsing them
  ## ahead of time:
  ##
  ##   1) Explicitly collapse the parent-child edges, check
  ##   that gap is zero.
  ##
  ##   2) Check other edges for a sample of vertices.  For
  ##   all negative gaps if not bfs, gather those edges and
  ##   re-check.

  pc_edge = ijw(parent_child_edge_list, :);
  [PC_i, PC_j, PC_w] = find (sparse (pc_edge(:, 1),
				     pc_edge(:, 2),
				     pc_edge(:, 3), N, N));
  if is_bfs,
     gap = 1 + d(PC_i) - d(PC_j);
  else
     gap = PC_w + d(PC_i) - d(PC_j);
  endif

  if any (gap != 0),
     ## Constraints not exactly satisfied along tree edges.
     out = -4;
     return;
  endif

  ## Determine which vertices to check.
  check_i = sample_roots (N, 2*SCALE, prngidx);
  to_check = ismember (pc_edge(:, 1), check_i) | \
    ismember (pc_edge(:, 2), check_i);
  to_check = pc_edge (to_check, :);

  [ijw_i, ijw_j, ijw_w] = find (sparse (to_check(:, 1),
					to_check(:, 2),
					to_check(:, 3),
					N, N));
  ## Note: This checks every edge adjacent to the
  ## vertices.  Only needs to check up to the edge
  ## factor.

  if is_bfs,
    gap = 1 + d(ijw_i) - d(ijw_j);
    if any (abs (gap) > 1),
      ## Some edge crosses two levels down the
      ## tree, cannot be from a BFS.
      out = -6;
      return;
    endif
  else
    gap = ijw_w + d(ijw_i) - d(ijw_j);
  endif

  if any (gap < 0),
    ## Dual constraint violated.
    out = -5;
    return;
  endif
endfunction
