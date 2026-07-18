---
name: jj-flow
description: jujutsu (jj) command cheatsheet for this colocated jj+git Bioconductor repo. Use when doing version-control work here — starting a change, committing, pushing a branch, opening a PR against devel, syncing, or when a git command is the reflex and the jj equivalent is wanted.
---

# jj workflow for GEOquery (colocated jj + git)

This repo has both `.jj/` and `.git/`. Every git command still works and the
remote only ever sees plain git branches, so the Bioconductor workflow (branches
pushed to `origin`, PRs merged into `devel`) is unchanged. `jj` just replaces the
local mechanics. `trunk()` is aliased to `devel@origin`; jj **bookmarks** are git
branches.

## Mental model

- The working copy **is a commit** (`@`). No staging area, no `git add`. Edits are
  auto-snapshotted into `@` on every `jj` command.
- You don't need a branch to start. Work anonymously; name a **bookmark** only
  when it's time to push.

## git → jj cheatsheet

| Task | git | jj |
|------|-----|-----|
| Start a change off devel | `git switch -c fix/x devel` | `jj new devel` |
| See state | `git status` | `jj st` |
| History | `git log --oneline` | `jj log` |
| Set/commit message | `git commit -m` | `jj describe -m "..."` (message on `@`) |
| Close change, start next | — | `jj commit -m "..."` |
| Amend current change | `git commit --amend` | just edit files (auto-snapshot) |
| Push a new branch (create+track) | `git push -u origin fix/x` | `jj git push --named fix/x=@` |
| Push updates | `git push` | `jj git push` |
| Fetch | `git fetch` | `jj git fetch` |
| Update onto latest devel | `git rebase devel` | `jj rebase -o devel` |
| Point an existing bookmark at a rev | `git branch -f fix/x @` | `jj bookmark set fix/x -r @` |
| Discard working changes | `git restore .` | `jj abandon` (drops `@`) or `jj restore` |
| Undo the last jj operation | `git reflog` + reset | `jj undo` |
| Rewind repo to an earlier state | — | `jj op restore <id>` |

## Typical PR flow

```sh
jj git fetch                                   # get latest devel@origin
jj new devel                                   # fresh change on top of trunk
# ...edit R/, tests, NEWS.md, bump DESCRIPTION...
jj describe -m "fix: parse GSExxxxx characteristics (#123)"
jj git push --named fix/issue-123=@               # creates, tracks, and pushes the bookmark
# open PR against devel on GitHub; merge when checks are green
```

Later updates to the same PR: edit, then `jj git push` (the bookmark auto-follows
`@` if you keep it there, or `jj bookmark set fix/issue-123 -r @` first).

Multi-commit PR: `jj commit -m "..."` closes the current change and opens a new
one on top; point the bookmark at the tip (`jj bookmark set fix/issue-123 -r @-`)
before pushing.

## Operation log (the undo button)

jj records every operation (commit, rebase, push, abandon…) in an **operation
log**, so almost nothing is unrecoverable:

- `jj op log` — list operations, newest first, each with an id.
- `jj undo` — revert the most recent operation.
- `jj op restore <id>` — rewind the *entire repo* to how it looked at that op
  (working copy, bookmarks, and all). This is the escape hatch after a bad
  rebase or an accidental `jj abandon`.

## Notes for this repo

- Conventional-commit prefixes still apply (`fix:` / `feat:` / `chore:` / `docs:`
  / `refactor:`) — see CONTRIBUTING.md.
- Commit trailer convention for AI-assisted commits still holds; add it in the
  `jj describe` message body.
- `devel` is branch-protected on GitHub — you can't push straight to it; always
  go through a bookmark + PR.
- Identity is set repo-locally in `.jj/` (not committed). If commits show an empty
  author, run `jj metaedit --update-author` after
  `jj config set --repo user.name/user.email`.
- Escape hatch: colocated means `git switch`, `git rebase`, `gh pr create`, etc.
  all still work on the same objects if you'd rather drop to git for a step.
