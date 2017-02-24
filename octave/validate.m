function out = validate (parent, ijw, search_key, d, is_sssp)
  %% Validate the results of BFS or SSSP.

  %% Default: no error.
  out = 1;

  %% Adjust from zero labels.
  parent = parent + 1;
  search_key = search_key + 1;
  ijw = ijw + 1;

  %% Remove self-loops.
  ijw(:,(ijw(1, :) == ijw(2, :))') = [];
  
  %% root must be the parent of itself.
  if parent (search_key) != search_key,
    out = 0;
    return;
  end

  N = max (max (ijw(1:2,:)));
  slice = find (parent > 0);

  %% Compute levels and check for cycles.
  level = zeros (size (parent));
  level (slice) = 1;
  P = parent (slice);
  mask = P != search_key;
  k = 0;
  while any (mask),
    level(slice(mask)) = level(slice(mask)) + 1;
    P = parent (P);
    mask = P != search_key;
    k = k + 1;
    if k > N,
      %% There must be a cycle in the tree.
      out = -3;
      return;
    end
  end

  %% Check that there are no edges with only one end in the tree.
  %% This also checks the component condition.
  lij = level (ijw(1:2,:));
  neither_in = lij(1,:) == 0 & lij(2,:) == 0;
  both_in = lij(1,:) > 0 & lij(2,:) > 0;
  if any (not (neither_in | both_in)),
    out = -4;
    return
  end

  %% Validate the distances/levels.
  respects_tree_level = true(1,size(ijw, 2));
  if !is_sssp
    respects_tree_level = abs (lij(1,:) - lij(2,:)) <= 1;
  else
    respects_tree_level = abs (d(ijw(1,:)) - d(ijw(2,:)))' <= ijw(3,:);
  end
  if any (not (neither_in | respects_tree_level))
    out = -5;
    return
  end
