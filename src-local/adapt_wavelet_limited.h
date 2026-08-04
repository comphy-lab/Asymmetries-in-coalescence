/**
# Position-dependent wavelet refinement ceiling

`MLFun(x,y,z)` supplies the maximum level for each cell. This is Pairetti's
regional-ceiling modification applied to the current Basilisk
`adapt_wavelet()` implementation. It lets the drop-injection drill retain full
resolution only in the end-pinchoff region while regularising unrelated
small-bubble-side singularities.
*/

trace
astats adapt_wavelet_limited (scalar * slist, double * max,
                              int (*MLFun)(double,double,double),
                              int minlevel = 1, scalar * list = all)
{
  scalar * ilist = list;
  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary (list);
    restriction (slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy ({cm, fm});
      for (scalar s in all)
        list = list_add (list, s);
    }
    boundary (list);
    scalar * listr = list_concat (slist, {cm});
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in list)
    listc = list_add_depend (listc, s);
  if (minlevel < 1)
    minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    int cellMAX = MLFun (x, y, z);
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
        if (cell.flags & too_coarse) {
          cell.flags &= ~too_coarse;
          refine_cell (point, listc, refined, &tree->refined);
          st.nf++;
        }
        continue;
      }
      else {
        if (cell.flags & refined) {
          cell.flags &= ~too_coarse;
          continue;
        }
        bool local = is_local(cell);
        if (!local)
          foreach_child()
            if (is_local(cell)) {
              local = true; break;
            }
        if (local) {
          int i = 0;
          static const int just_fine = 1 << (user + 3);
          for (scalar s in slist) {
            double emax = max[i++], sc[(1 << dimension)*s.block];
            double * b = sc;
            foreach_child()
              foreach_blockf(s)
                *b++ = s[];
            s.prolongation (point, s);
            b = sc;
            foreach_child()
              foreach_blockf(s) {
                double e = fabs(*b - s[]);
                if (e > emax && level < cellMAX) {
                  cell.flags &= ~too_fine;
                  cell.flags |= too_coarse;
                }
                else if ((e <= emax/1.5 || level > cellMAX) &&
                         !(cell.flags & (too_coarse|just_fine))) {
                  if (level >= minlevel)
                    cell.flags |= too_fine;
                }
                else if (!(cell.flags & too_coarse)) {
                  cell.flags &= ~too_fine;
                  cell.flags |= just_fine;
                }
                s[] = *b++;
              }
          }
          foreach_child() {
            cell.flags &= ~just_fine;
            if (!is_leaf(cell)) {
              cell.flags &= ~too_coarse;
              if (level >= cellMAX)
                cell.flags |= too_fine;
            }
            else if (!is_active(cell))
              cell.flags &= ~too_coarse;
          }
        }
      }
    }
  }
  mpi_boundary_refine (listc);

  for (int l = depth(); l >= 0; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
        if (level == l) {
          if (!is_leaf(cell)) {
            if (cell.flags & refined)
              cell.flags &= ~(refined|too_fine);
            else if (cell.flags & too_fine) {
              if (is_local(cell) && coarsen_cell (point, listc))
                st.nc++;
              cell.flags &= ~too_fine;
            }
          }
          if (cell.flags & too_fine)
            cell.flags &= ~too_fine;
          else if (level > 0 && (aparent(0).flags & too_fine))
            aparent(0).flags &= ~too_fine;
          continue;
        }
        else if (is_leaf(cell))
          continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);
  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);
  if (list != ilist)
    free (list);
  return st;
}
