#include "internal.hpp"

namespace CaDiCaL {

/*------------------------------------------------------------------------*/

// Failed literal probing uses it's own propagation and assignment
// functions.  It further provides on-the-fly generation of hyper binary
// resolvents but only probes on roots of the binary implication graph.  The
// search for failed literals is limited, but untried roots are kept until
// the next time 'probe' is called.  Left over probes from the last attempt
// and new probes are tried until the limit is hit or all are tried.

/*------------------------------------------------------------------------*/

bool Internal::probing () {
  if (!opts.probe) return false;
  if (!preprocessing && !opts.inprocessing) return false;
  if (preprocessing) assert (lim.preprocessing);
  if (stats.probingphases &&
      last.probe.reductions == stats.reductions) return false;
  return lim.probe <= stats.conflicts;
}

/*------------------------------------------------------------------------*/

inline int Internal::get_parent_reason_literal (int lit) {
  const int idx = vidx (lit);
  int res = parents[idx];
  if (lit < 0) res = -res;
  return res;
}

inline void Internal::set_parent_reason_literal (int lit, int reason) {
  const int idx = vidx (lit);
  if (lit < 0) reason = -reason;
  parents[idx] = reason;
}

/*------------------------------------------------------------------------*/

// On-the-fly (dynamic) hyper binary resolution on decision level one can
// make use of the fact that the implication graph is actually a tree.

// Compute a dominator of two literals in the binary implication tree.

int Internal::probe_dominator (int a, int b) {
  require_mode (PROBE);
  int l = a, k = b;
  Var * u = &var (l), * v = &var (k);
  assert (val (l) > 0), assert (val (k) > 0);
  assert (u->level == 1), assert (v->level == 1);
  while (l != k) {
    if (u->trail > v->trail) swap (l, k), swap (u, v);
    if (!get_parent_reason_literal (l)) return l;
    int parent = get_parent_reason_literal (k);
    assert  (parent), assert (val (parent) > 0);
    v = &var (k = parent);
    assert (v->level == 1);
  }
  LOG ("dominator %d of %d and %d", l, a, b);
  assert (val (l) > 0);
  return l;
}

// The idea of dynamic on-the-fly hyper-binary resolution came up in the
// PrecoSAT solver, where it originally was used on all decision levels.

// It turned out, that most of the hyper-binary resolvents were generated
// during probing on decision level one anyhow.  Thus this version is
// specialized to decision level one, where actually all long (non-binary)
// forcing clauses can be resolved to become binary.  So if we find a clause
// which would force a new assignment at decision level one during probing
// we resolve it (the 'reason' argument) to obtain a hyper binary resolvent.
// It consists of the still unassigned literal (the new unit) and the
// negation of the unique closest dominator of the negation of all (false)
// literals in the clause (which has to exist on decision level one).

// There are two special cases which should be mentioned:
//
//   (A) The reason is already a binary clause in a certain sense, since all
//   its unwatched literals are root level fixed to false.  In this
//   situation it would be better to shrink the clause immediately instead
//   of adding a new clause consisting only of the watched literals.
//   However, this would happen during the next garbage collection anyhow.
//
//   (B) The resolvent subsumes the original reason clause. This is
//   equivalent to the property that the negated dominator is contained in
//   the original reason.  Again one could in principle shrink the clause.
//
// Note that (A) is actually subsumed by (B).  The possible optimization to
// shrink the clause on-the-fly is difficult (need to update 'blit' and
// 'binary' of the other watch at least) and also not really that important.
// For (B) we simply add the new binary resolvent and mark the old subsumed
// clause as garbage instead.  And since in the situation of (A) the
// shrinking will be performed at the next  garbage collection anyhow, we
// do not change clauses in (A).

// The hyper binary resolvent clause is redundant unless it subsumes the
// original reason and that one is irredundant.

// If the option 'opts.probehbr' is 'false', we actually do not add the new
// hyper binary resolvent, but simply pretend we would have added it and
// still return the dominator as new reason / parent for the new unit.

// Finally note that adding clauses changes the watches of the propagated
// literal and thus we can not use standard iterators during probing but
// need to fall back to indices.  One watch for the hyper binary resolvent
// clause is added at the end of the currently propagated watches, but its
// watch is a binary watch and will be skipped during propagating long
// clauses anyhow.

inline int Internal::hyper_binary_resolve (Clause * reason) {
  require_mode (PROBE);
  assert (level == 1);
  assert (reason->size > 2);
  const const_literal_iterator end = reason->end ();
  const int * lits = reason->literals;
  const_literal_iterator k;
#ifndef NDEBUG
  // First literal unassigned, all others false.
  assert (!val (lits[0]));
  for (k = lits + 1; k != end; k++) assert (val (*k) < 0);
  assert (var (lits[1]).level == 1);
#endif
  LOG (reason, "hyper binary resolving");
  stats.hbrs++;
  stats.hbrsizes += reason->size;
  const int lit = lits[1];
  int dom = -lit, non_root_level_literals = 0;
  for (k = lits + 2; k != end; k++) {
    const int other = -*k;
    assert (val (other) > 0);
    if (!var (other).level) continue;
    dom = probe_dominator (dom, other);
    non_root_level_literals++;
  }
  if (non_root_level_literals && opts.probehbr) { // !(A)
    bool contained = false;
    for (k = lits + 1; !contained && k != end; k++)
      contained = (*k == -dom);
    const bool red = !contained || reason->redundant;
    if (red) stats.hbreds++;
    LOG ("new %s hyper binary resolvent %d %d",
      (red ? "redundant" : "irredundant"), -dom, lits[0]);
    assert (clause.empty ());
    clause.push_back (-dom);
    clause.push_back (lits[0]);
    Clause * c = new_hyper_binary_resolved_clause (red, 2);
    if (red) c->hyper = true;
    clause.clear ();
    if (contained) {
      stats.hbrsubs++;
      LOG (reason, "subsumed original");
      mark_garbage (reason);
    }
  }
  return dom;
}

/*------------------------------------------------------------------------*/

// The following functions 'probe_assign' and 'probe_propagate' are used for
// propagating during failed literal probing in simplification mode, as
// replacement of the generic propagation routine 'propagate' and
// 'search_assign'.

// The code is mostly copied from 'propagate.cpp' and specialized.  We only
// comment on the differences.  More explanations are in 'propagate.cpp'.

inline void Internal::probe_assign (int lit, int parent) {
  require_mode (PROBE);
  int idx = vidx (lit);
  assert (!vals[idx]);
  assert (!flags (idx).eliminated () || !parent);
  assert (!parent || val (parent) > 0);
  Var & v = var (idx);
  v.level = level;
  v.trail = (int) trail.size ();
  set_parent_reason_literal (lit, parent);
  if (!level) learn_unit_clause (lit);
  else assert (level == 1);
  const signed char tmp = sign (lit);
  vals[idx] = tmp;
  vals[-idx] = -tmp;
  assert (val (lit) > 0);
  assert (val (-lit) < 0);
  trail.push_back (lit);

  // Do not save the current phase during inprocessing but remember the
  // number of units on the trail of the last time this literal was
  // assigned.  This allows us to avoid some redundant failed literal
  // probing attempts.  Search for 'propfixed' in 'probe.cpp' for details.
  //
  if (level) propfixed (lit) = stats.all.fixed;

      if (parent) LOG ("probe assign %d parent %d", lit, parent);
  else if (level) LOG ("probe assign %d probe", lit);
  else            LOG ("probe assign %d negated failed literal UIP", lit);
}

void Internal::probe_assign_decision (int lit) {
  require_mode (PROBE);
  assert (!level);
  assert (propagated == trail.size ());
  level++;
  control.push_back (Level (lit, trail.size ()));
  probe_assign (lit, 0);
}

void Internal::probe_assign_unit (int lit) {
  require_mode (PROBE);
  assert (!level);
  assert (active (lit));
  probe_assign (lit, 0);
}

/*------------------------------------------------------------------------*/

// This is essentially the same as 'propagate' except that we prioritize and
// always propagate binary clauses first (see our CPAIOR'13 paper on tree
// based look ahead), then immediately stop at a conflict and of course use
// 'probe_assign' instead of 'search_assign'.  The binary propagation part
// is factored out too.  If a new unit on decision level one is found we
// perform hyper binary resolution and thus actually build an implication
// tree instead of a DAG.  Statistics counters are also different.

inline void Internal::probe_propagate2 () {
  require_mode (PROBE);
  while (propagated2 != trail.size ()) {
    const int lit = -trail[propagated2++];
    LOG ("probe propagating %d over binary clauses", -lit);
    Watches & ws = watches (lit);
    for (const auto & w : ws) {
      if (!w.binary ()) continue;
      const signed char b = val (w.blit);
      if (b > 0) continue;
      if (b < 0) conflict = w.clause;                   // but continue
      else probe_assign (w.blit, -lit);
    }
  }
}

bool Internal::probe_propagate () {
  require_mode (PROBE);
  assert (!unsat);
  START (propagate);
  int64_t before = propagated2 = propagated;
  while (!conflict) {
    if (propagated2 != trail.size ()) probe_propagate2 ();
    else if (propagated != trail.size ()) {
      const int lit = -trail[propagated++];
      LOG ("probe propagating %d over large clauses", -lit);
      Watches & ws = watches (lit);
      size_t i = 0, j = 0;
      while (i != ws.size ()) {
        const Watch w = ws[j++] = ws[i++];
        if (w.binary ()) continue;
        const signed char b = val (w.blit);
        if (b > 0) continue;
        if (w.clause->garbage) continue;
        const literal_iterator lits = w.clause->begin ();
        const int other = lits[0]^lits[1]^lit;
        //lits[0] = other, lits[1] = lit;
        const signed char u = val (other);
        if (u > 0) ws[j-1].blit = other;
        else {
          const int size = w.clause->size;
          const const_literal_iterator end = lits + size;
          const literal_iterator middle = lits + w.clause->pos;
          literal_iterator k = middle;
          int r = 0;
          signed char v = -1;
          while (k != end && (v = val (r = *k)) < 0)
            k++;
          if (v < 0) {
            k = lits + 2;
            assert (w.clause->pos <= size);
            while (k != middle && (v = val (r = *k)) < 0)
              k++;
          }
          w.clause->pos = k - lits;
          assert (lits + 2 <= k), assert (k <= w.clause->end ());
          if (v > 0) ws[j-1].blit = r;
          else if (!v) {
            LOG (w.clause, "unwatch %d in", r);
            *k = lit;
            lits[0] = other;
            lits[1] = r;
            watch_literal (r, lit, w.clause);
            j--;
          } else if (!u) {
            if (level == 1) {
              lits[0] = other, lits[1] = lit;
              int dom = hyper_binary_resolve (w.clause);
              probe_assign (other, dom);
            } else probe_assign_unit (other);
            probe_propagate2 ();
          } else conflict = w.clause;
        }
      }
      if (j != i) {
        while (i != ws.size ())
          ws[j++] = ws[i++];
        ws.resize (j);
      }
    } else break;
  }
  int64_t delta = propagated2 - before;
  stats.propagations.probe += delta;
  if (conflict) LOG (conflict, "conflict");
  STOP (propagate);
  return !conflict;
}

/*------------------------------------------------------------------------*/

// This a specialized instance of 'analyze'.

void Internal::failed_literal (int failed) {

  LOG ("analyzing failed literal probe %d", failed);
  stats.failed++;
  stats.probefailed++;

  assert (!unsat);
  assert (conflict);
  assert (level == 1);
  assert (analyzed.empty ());

  START (analyze);

  LOG (conflict, "analyzing failed literal conflict");

  int uip = 0;
  for (const auto & lit : *conflict) {
    const int other = -lit;
    if (!var (other).level) continue;
    uip = uip ? probe_dominator (uip, other) : other;
  }
  LOG ("found probing UIP %d", uip);
  assert (uip);

  vector<int> work;
  int parent = uip;
  while (parent != failed) {
    const int next = get_parent_reason_literal (parent);
    parent = next;
    assert (parent);
    work.push_back (parent);
  }

  backtrack ();
  clear_analyzed_literals ();
  conflict = 0;

  assert (!val (uip));
  probe_assign_unit (-uip);

  if (!probe_propagate ()) learn_empty_clause ();

  while (!unsat && !work.empty ()) {
    const int parent = work.back ();
    work.pop_back ();
    const signed char tmp = val (parent);
    if (tmp < 0) continue;
    if (tmp > 0) {
      LOG ("clashing failed parent %d", parent);
      learn_empty_clause ();
    } else {
      LOG ("found unassigned failed parent %d", parent);
      probe_assign_unit (-parent);
      if (!probe_propagate ()) learn_empty_clause ();
    }
  }
  erase_vector (work);

  STOP (analyze);

  assert (unsat || val (failed) < 0);
}

/*------------------------------------------------------------------------*/

bool Internal::is_binary_clause (Clause * c, int & a, int & b) {
  assert (!level);
  if (c->garbage) return false;
  int first = 0, second = 0;
  for (const auto & lit : *c) {
    const signed char tmp = val (lit);
    if (tmp > 0) return false;
    if (tmp < 0) continue;
    if (second) return false;
    if (first) second = lit;
    else first = lit;
  }
  if (!second) return false;
  a = first, b = second;
  return true;
}

// We probe on literals first, which occur more often negated and thus we
// sort the 'probes' stack in such a way that literals which occur negated
// less frequently come first.  Probes are taken from the back of the stack.

struct probe_negated_noccs_rank {
  Internal * internal;
  probe_negated_noccs_rank (Internal * i) : internal (i) { }
  size_t operator () (int a) const { return internal->noccs (-a); }
};

// Fill the 'probes' schedule.

void Internal::generate_probes () {

  assert (probes.empty ());

  // First determine all the literals which occur in binary clauses. It is
  // way faster to go over the clauses once, instead of walking the watch
  // lists for each literal.
  //
  init_noccs ();
  for (const auto & c : clauses) {
    int a, b;
    if (!is_binary_clause (c, a, b)) continue;
    noccs (a)++;
    noccs (b)++;
  }

  for (auto idx : vars) {

    // Then focus on roots of the binary implication graph, which are
    // literals occurring negatively in a binary clause, but not positively.
    // If neither 'idx' nor '-idx' is a root it makes less sense to probe
    // this variable.

    // This argument requires that equivalent literal substitution through
    // 'decompose' is performed, because otherwise there might be 'cyclic
    // roots' which are not tried, i.e., -1 2 0, 1 -2 0, 1 2 3 0, 1 2 -3 0.

    const bool have_pos_bin_occs = noccs (idx) > 0;
    const bool have_neg_bin_occs = noccs (-idx) > 0;

    if (have_pos_bin_occs == have_neg_bin_occs) continue;

    int probe = have_neg_bin_occs ? idx : -idx;

    // See the discussion where 'propfixed' is used below.
    //
    if (propfixed (probe) >= stats.all.fixed) continue;

    LOG ("scheduling probe %d negated occs %" PRId64 "", probe, noccs (-probe));
    probes.push_back (probe);
  }

  rsort (probes.begin (), probes.end (), probe_negated_noccs_rank (this));

  reset_noccs ();
  shrink_vector (probes);

  PHASE ("probe-round", stats.probingrounds,
    "scheduled %" PRId64 " literals %.0f%%",
    probes.size (), percent (probes.size (), 2u*max_var));
}

// Follow the ideas in 'generate_probes' but flush non root probes and
// reorder remaining probes.

void Internal::flush_probes () {

  assert (!probes.empty ());

  init_noccs ();
  for (const auto & c : clauses) {
    int a, b;
    if (!is_binary_clause (c, a, b)) continue;
    noccs (a)++;
    noccs (b)++;
  }

  const auto eop = probes.end ();
  auto j = probes.begin ();
  for (auto i = j; i != eop; i++) {
    int lit = *i;
    if (!active (lit)) continue;
    const bool have_pos_bin_occs = noccs (lit) > 0;
    const bool have_neg_bin_occs = noccs (-lit) > 0;
    if (have_pos_bin_occs == have_neg_bin_occs) continue;
    if (have_pos_bin_occs) lit = -lit;
    assert (!noccs (lit)), assert (noccs (-lit) > 0);
    if (propfixed (lit) >= stats.all.fixed) continue;
    LOG ("keeping probe %d negated occs %" PRId64 "", lit, noccs (-lit));
    *j++ = lit;
  }
  size_t remain = j - probes.begin ();
#ifndef QUIET
  size_t flushed = probes.size () - remain;
#endif
  probes.resize (remain);

  rsort (probes.begin (), probes.end (), probe_negated_noccs_rank (this));

  reset_noccs ();
  shrink_vector (probes);

  PHASE ("probe-round", stats.probingrounds,
    "flushed %zd literals %.0f%% remaining %zd",
    flushed, percent (flushed, remain + flushed), remain);
}

int Internal::next_probe () {

  int generated = 0;

  for (;;) {

    if (probes.empty ()) {
      if (generated++) return 0;
      generate_probes ();
    }

    while (!probes.empty ()) {

      int probe = probes.back ();
      probes.pop_back ();

      // Eliminated or assigned.
      //
      if (!active (probe)) continue;

      // There is now new unit since the last time we propagated this probe,
      // thus we propagated it before without obtaining a conflict and
      // nothing changed since then.  Thus there is no need to propagate it
      // again.  This observation was independently made by Partik Simons
      // et.al. in the context of implementing 'smodels' (see for instance
      // Alg. 4 in his JAIR article from 2002) and it has also been
      // contributed to the thesis work of Yacine Boufkhad.
      //
      if (propfixed (probe) >= stats.all.fixed) continue;

      return probe;
    }
  }
}

bool Internal::probe_round () {

  if (unsat) return false;
  if (terminated_asynchronously ()) return false;

  START_SIMPLIFIER (probe, PROBE);
  stats.probingrounds++;

  // Probing is limited in terms of non-probing propagations
  // 'stats.propagations'. We allow a certain percentage 'opts.probereleff'
  // (say %5) of probing propagations in each probing with a lower bound of
  // 'opts.probmineff'.
  //
  int64_t delta = stats.propagations.search;
  delta -= last.probe.propagations;
  delta *= 1e-3 * opts.probereleff;
  if (delta < opts.probemineff) delta = opts.probemineff;
  if (delta > opts.probemaxeff) delta = opts.probemaxeff;
  delta += 2l * active ();

  PHASE ("probe-round", stats.probingrounds,
    "probing limit of %" PRId64 " propagations ", delta);

  int64_t limit = stats.propagations.probe + delta;

  int old_failed = stats.failed;
#ifndef QUIET
  int64_t old_probed = stats.probed;
#endif
  int64_t old_hbrs = stats.hbrs;

  if (!probes.empty ()) flush_probes ();

  // We reset 'propfixed' since there was at least another conflict thus
  // a new learned clause, which might produce new propagations (and hyper
  // binary resolvents).  During 'generate_probes' we keep the old value.
  //
  for (auto idx : vars)
    propfixed (idx) = propfixed (-idx) = -1;

  assert (unsat || propagated == trail.size ());
  propagated = propagated2 = trail.size ();

  int probe;
  while (!unsat &&
         !terminated_asynchronously () &&
         stats.propagations.probe < limit &&
         (probe = next_probe ())) {
    stats.probed++;
    LOG ("probing %d", probe);
    probe_assign_decision (probe);
    if (probe_propagate ()) backtrack ();
    else failed_literal (probe);
  }

  if (unsat) LOG ("probing derived empty clause");
  else if (propagated < trail.size ()) {
    LOG ("probing produced %" PRId64 " units", trail.size () - propagated);
    if (!propagate ()) {
      LOG ("propagating units after probing results in empty clause");
      learn_empty_clause ();
    } else sort_watches ();
  }

  int failed = stats.failed - old_failed;
#ifndef QUIET
  int64_t probed = stats.probed - old_probed;
#endif
  int64_t hbrs = stats.hbrs - old_hbrs;

  PHASE ("probe-round", stats.probingrounds,
    "probed %" PRId64 " and found %d failed literals", probed, failed);

  if (hbrs)
    PHASE ("probe-round", stats.probingrounds,
      "found %" PRId64 " hyper binary resolvents", hbrs);

  STOP_SIMPLIFIER (probe, PROBE);

  report ('p', !opts.reportall && !(unsat + failed + hbrs));

  return !unsat && failed;
}

/*------------------------------------------------------------------------*/

void CaDiCaL::Internal::probe (bool update_limits) {

  if (unsat) return;
  if (level) backtrack ();
  if (!propagate ()) { learn_empty_clause (); return; }

  stats.probingphases++;

  const int before = active ();

  // We trigger equivalent literal substitution (ELS) before ...
  //
  decompose ();

  if (ternary ())       // If we derived a binary clause
    decompose ();       // then start another round of ELS.

  // Remove duplicated binary clauses and perform in essence hyper unary
  // resolution, i.e., derive the unit '2' from '1 2' and '-1 2'.
  //
  mark_duplicated_binary_clauses_as_garbage ();

  for (int round = 1; round <= opts.proberounds; round++)
    if (!probe_round ())
      break;

  decompose ();         // ... and (ELS) afterwards.

  last.probe.propagations = stats.propagations.search;

  if (!update_limits) return;

  const int after = active ();
  const int removed = before - after;
  assert (removed >= 0);

  if (removed) {
    stats.probesuccess++;
    PHASE ("probe-phase", stats.probingphases,
      "successfully removed %d active variables %.0f%%",
      removed, percent (removed, before));
  } else
    PHASE ("probe-phase", stats.probingphases,
      "could not remove any active variable");

  const int64_t delta = opts.probeint * (stats.probingphases + 1);
  lim.probe = stats.conflicts + delta;

  PHASE ("probe-phase", stats.probingphases,
    "new limit at %" PRId64 " conflicts after %" PRId64 " conflicts",
    lim.probe, delta);

  last.probe.reductions = stats.reductions;
}

}
