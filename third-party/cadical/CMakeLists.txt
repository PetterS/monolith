add_definitions("-DNBUILD")
if (MSVC)
    add_definitions("-DNUNLOCKED")
endif ()

add_library(cadical
cadical/analyze.cpp    cadical/decompose.cpp    cadical/minimize.cpp   cadical/signal.cpp
cadical/arena.cpp      cadical/deduplicate.cpp  cadical/occs.cpp       cadical/solution.cpp
cadical/assume.cpp     cadical/elim.cpp         cadical/options.cpp    cadical/solver.cpp
cadical/averages.cpp   cadical/ema.cpp          cadical/parse.cpp      cadical/stats.cpp
cadical/backtrack.cpp  cadical/extend.cpp       cadical/phases.cpp     cadical/subsume.cpp
cadical/backward.cpp   cadical/external.cpp     cadical/probe.cpp      cadical/terminal.cpp
cadical/bins.cpp       cadical/file.cpp         cadical/profile.cpp    cadical/ternary.cpp
cadical/block.cpp      cadical/flags.cpp        cadical/proof.cpp      cadical/tracer.cpp
cadical/ccadical.cpp   cadical/format.cpp       cadical/propagate.cpp  cadical/transred.cpp
cadical/checker.cpp    cadical/gates.cpp        cadical/queue.cpp      cadical/util.cpp
cadical/clause.cpp     cadical/instantiate.cpp  cadical/random.cpp     cadical/var.cpp
cadical/collect.cpp    cadical/internal.cpp     cadical/reduce.cpp     cadical/version.cpp
cadical/compact.cpp    cadical/ipasir.cpp       cadical/rephase.cpp    cadical/vivify.cpp
cadical/condition.cpp  cadical/limit.cpp        cadical/report.cpp     cadical/walk.cpp
cadical/config.cpp     cadical/logging.cpp      cadical/resources.cpp  cadical/watch.cpp
cadical/contract.cpp   cadical/lookahead.cpp    cadical/restart.cpp
cadical/cover.cpp      cadical/lucky.cpp        cadical/restore.cpp
cadical/decide.cpp     cadical/message.cpp      cadical/score.cpp)

