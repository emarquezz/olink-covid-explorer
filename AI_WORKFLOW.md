# AI_WORKFLOW — How this project was built with AI assistance

This document records how AI assistance was used in the Olink COVID-19
Explorer project: the working protocol, and — more importantly — an
honest log of what failed and what each failure taught. The failure log
is the point. Anyone can claim a smooth AI-assisted project; the
verifiable value is in the loop that caught the failures.

## The core loop

Every unit of work followed one loop, without exception:

  prompt → code → run it → paste console output back → verify → commit

- **The AI writes code; the human runs it.** The AI never executed
anything. Every artifact, error message, and console line quoted in
this repo was pasted back by the human from a real run.
- **The console output is the evidence chain.** Nothing entered the
README or report without a pasted console line behind it. Numbers
like "436 proteins × 307 samples, 0 dropped, 275 imputed" trace to
run logs, not to AI assertions.
- **Regression checks every step.** After every change, the previously
working features were re-verified ("the other tabs still work").
This repeatedly caught failures the change had no business causing
(see Anecdote #4).
  - **One change per commit**, with summary + description. The commit
  history doubles as the project's decision log, including honest
  notes about fixes and iterations.

## Role split

| AI was good at | AI was bad at (and the loop caught it) |
|---|---|
| Drafting complete pipelines with consistent structure | Full-file consolidation — corrupted untouched code (Anecdote #4) |
| Explaining error messages and proposing hypotheses | Confident first hypotheses that were wrong (Anecdote #3) |
| Refactoring toward the artifact contract (`result` + `params`) | Output-size limits — silent truncation of long deliveries (Anecdote #6) |
| Turning verified console facts into documentation prose | Judging legibility of rendered UI it could not see (Anecdote #5) |

The pattern: AI is trusted as a **draftsman and analyst**, never as the
**executor or verifier**. Verification always ran through the human's
  console and eyes.

  ---

    # Failure log

    ## Anecdote #1 — The error that pointed at the wrong place

    The first `harmonize()` run failed on join columns, reporting:

    Optional column absent: WHO_Severity_Peak

  — for the **plasma** subcohort, which demonstrably *has* that column.
  The message was both true and misleading: an earlier step had already
  dropped the column by the time the check ran.

  **Root cause:** two chained `standardize()` calls. The first call's
`select()` silently dropped the optional column; the second call then
overwrote `meta` with *only* that column — each step was individually
plausible, the composition was broken.

**Lessons:**
1. Error location ≠ error cause. Messages must be read in execution
   order, as a chain — not in isolation.
2. When a message contradicts the data ("column absent" for a column
   that exists), suspect the pipeline's intermediate state, not the
  data.
  3. The fix was structural, not patched: one `standardize()` call with
  an `optional` argument, instead of composing two lossy calls.

  ## Anecdote #2 — Two errors, one per R dialect, in one function

  `run_de()` failed twice, in sequence:

    1. `all_of(covariates) must be size 133577 or 1, not 2` — inside a
  `distinct()` call.
  2. `Age(40,60] … non-valid names` — binned factor levels leaking into
  `model.matrix()` names, which `makeContrasts()` then tried to parse
  as R code.

  **Root causes:** two R-dialect collisions, unrelated to each other:
    (a) `all_of()` behaves as a tidyselect helper inside `select()` but as
  a data-mask *expression* inside `distinct()` — same function, opposite
  semantics, depending on context; (b) `makeContrasts()` parses its
  argument as R source, so any factor level that isn't a valid R name
(comma-containing Age bins like `Age(40,60]`) explodes.

**Fix:** `select(any_of(...))` instead of `all_of()` where optional;
`distinct()` without the masking trap; and — decisively — **the
contrast vector built by hand** as a named numeric vector, instead of
parsing a contrast string through the R interpreter.

**Lessons:**
1. Tidyselect helpers are context-dependent — `all_of()` is not
   one function with one semantics. When identical syntax works in
   one verb and fails in another, suspect the evaluation context,
   not the syntax.
2. Pass *values* (a numeric contrast vector), not *code* (a string
   to be parsed), across package boundaries — it eliminates a whole
   failure class at the dplyr ↔ limma seam.

## Anecdote #3 — The wrong first hypothesis, and the minimal repro

The volcano tab crashed with:

    is.character(txt) is not TRUE

The first hypothesis — a ggplotly/ggplot2 incompatibility — was
**disproved by a minimal repro**: the smallest volcano render worked
fine outside the app. The traceback was then read properly, and the
failing frame was not in the plotting at all: it was the `validate()`
call inside `de_results()`.

**Root cause:** `plotly::validate` **masks** `shiny::validate** under
`library()` load order. The app was calling the wrong `validate()` —
plotly's, which checks `txt` — with a `need()` condition as argument.

  **Fix:** every validation call namespaced as
  `shiny::validate(shiny::need(...))`. The fix survives in the shipped
  code as a comment block documenting exactly this masking hazard.

  **Lessons:**
    1. Read the full traceback **before** theorizing; the failing frame
  names the culprit, and it wasn't where the symptom looked like it
   lived.
2. A minimal repro is the fastest way to *disprove* a hypothesis —
   its value is negative results.
3. Namespace any function name that exists in both shiny and plotly
   (`validate`, `layout`, …) — masking is silent until it isn't.

  ## Anecdote #4 — Consolidation drift (the volcano regression)

  Step 12 added a Heatmap tab — new code, merged into an app that worked.
  The step delivered a full-file rewrite. The app then crashed on the
  volcano:

    Error in data.frame: arguments imply differing number of rows: 3, 1, 0

  — inside a block that had not logically changed since the volcano was
  first written. The error message was a perfect column census: three
  vectors of length 3, one of length 1, one of length 0 — the corrupted
  `data.frame()` had a scalar `yend` and a stray NULL-named column where
  one clean line used to be.

  **Root cause:** merge-by-rewrite. Assembling the new file corrupted an
  untouched section; nothing about the volcano step was supposed to
  change, which made the bug instantly localizable: diff the new file
  against the last working version — the delta is the bug.

  **Fix:** the correct block restored from the last working version. All
  subsequent merge steps switched to surgical patches: paste-in blocks
  with stable output IDs, never full-file rewrites.

  **Lessons:**
    1. Full-file rewrites risk regressions in untouched sections. Diff
  against the last known-good file before running.
  2. The loop's regression check — "the other tabs still work" — is what
   catches this failure class. It is not bureaucracy.
3. Read error messages literally: "3, 1, 0" named exactly which
   vectors disagreed, in order.

## Anecdote #5 — "It rendered" is not "it's readable"

The QC tab first shipped as a 2x2 `layout_columns` grid — four cards,
because the metrics came in fours. The grid sizes each row to its
tallest content; the two plotly widgets resized to their squeezed
containers and rendered as unreadable slivers. No error, no warning,
no failed check: the app "worked" and the figures were invisible.

**Root cause:** layout chosen by widget count, not by content.
Searchable 10-column DTs and bar charts with 20 rotated labels both
want full width.

**Fix:** the tab restructured into full-width nested `navset_card_tab`
sub-tabs — a UI-only change, zero server edits, possible only because
the surgical-patch discipline had kept every output ID stable.

**Lessons:**
1. Content chooses layout.
2. Visual QA belongs in the loop: "it works" means "legible at real
  window sizes." A figure that renders illegibly is a bug.
3. Record the iteration honestly. The commit message states that the
   first layout squeezed the widgets and why sub-tabs replaced it —
   reviewers read that as judgment, not failure.

## Anecdote #6 — Chunked delivery as a protocol, not an improvisation

`report.qmd` could not be delivered in one message. The output
truncated silently, mid-file, repeatedly across multiple requests —
each time the file simply stopped, with no error and no warning that
anything was missing.

**Root cause:** hard output-size limits on the delivery channel. The
complete file exceeded them; truncation was silent, so "did the file
  end where it should?" became a legitimate check question.

**Fix:** chunking made a declared protocol rather than a rescue
attempt: explicit part boundaries (1, 2a, 2b), an assembly order
("append directly under the concordance section"), a "don't render
until complete" rule, and verification of the joined file end-to-end.

This is the symmetric twin of Anecdote #4. Both failures are "too
much in one unit" — the full-file rewrite corrupted one line; the
full-file delivery cut the file. Both were fixed the same way:
smaller units with explicit boundaries, then verification at the join.

**Lessons:**
1. When a channel has hard size limits, make chunking a first-class
   protocol: declared boundaries, assembly order, join verification.
2. The deliverable is the *joined* file — and joins need verification.
3. The last line of every part should be predictable; if it isn't
  what you expected, the part is incomplete.

  ## Retrospective — what the failure log shows

  | # | Failure | Root cause | Fix / lesson |
    |---|---|---|---|
    | 1 | `harmonize()` reported an optional column absent for plasma — which has it | two chained `standardize()` calls: first dropped the column, second overwrote `meta` | error location ≠ cause; structural fix — one call with an `optional` argument |
    | 2 | `run_de()`: `all_of()` size error, then invalid names from `Age(40,60]` | tidyselect vs data-mask semantics; `makeContrasts()` parses strings as R | `any_of()` where optional; hand-built numeric contrast vector — pass values, not code, across package boundaries |
    | 3 | Volcano: `is.character(txt) is not TRUE`; first hypothesis (ggplotly incompatibility) wrong | `plotly::validate` masks `shiny::validate` under library load order | read the full traceback before theorizing; minimal repro disproves; namespace colliding functions |
    | 4 | Volcano broke while adding an unrelated tab | full-file consolidation corrupted an untouched block | diff before running; prefer surgical patches |
    | 5 | QC figures invisible despite "rendering" | layout chosen by widget count, not content | content chooses layout; visual QA in the loop |
    | 6 | `report.qmd` truncated mid-delivery, repeatedly | hard output-size limit on the channel | declared part boundaries + join verification |

    The pattern across all six: every failure was caught by the same net —
  *run it, look at it, compare against the last known-good state.* None
  was caught by hoping. And no failure was hidden: each one is committed
  with its fix and its lesson, so the history and this document agree.

  ## The workflow, distilled

  1. **One loop, everywhere:** prompt → code → run → paste console output
  → verify → commit. Nothing is trusted that wasn't run; the console
   output *is* the evidence chain.
2. **A verified-facts ledger:** every number in the README and report
   traces to a pasted console line from a real run.
3. **Artifacts describe themselves:** `list(result, params)` — apps and
   reports read `params` and adapt; no hardcoded subcohort branching.
4. **Surgical patches over rewrites:** merge steps add blocks with
   stable IDs; full-file rewrites get diffed first.
5. **Chunked delivery with explicit joins:** large artifacts arrive as
   declared parts; the joined file is verified end-to-end.
6. **Visual QA is part of the loop:** a figure that renders illegibly
   is a bug, not a style preference.
7. **Honest history in commits:** iterations and fixes are recorded in
   commit descriptions — `git log --oneline` tells the same story this
   document does.
