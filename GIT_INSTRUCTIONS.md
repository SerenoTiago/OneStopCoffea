# Git instructions for the hadronic SUSY analysis

This repository uses two GitHub remotes and two branches with distinct roles:

- `origin`: Tiago's fork, `SerenoTiago/OneStopCoffea`.
- `upstream`: the original collaboration repository, `UMN-CMS/OneStopCoffea`.
- `master`: a clean local copy of `upstream/master`. Do not develop here.
- `analysis/hadronic-susy`: the working branch for the hadronic SUSY analysis.

Pushing to `upstream` is disabled locally as a safety measure. Analysis changes
should only be pushed to the fork through `origin`.

## Normal working routine

Start each work session from the analysis branch:

```bash
cd /uscms/home/tsereno/nobackup/OneStopCoffea
git switch analysis/hadronic-susy
git status
```

After completing a coherent change, inspect it before committing:

```bash
git status
git diff
```

Stage only the files that belong to that change:

```bash
git add path/to/file1 path/to/file2
git diff --staged
```

Then commit and push:

```bash
git commit -m "Briefly describe the completed change"
git push
```

Prefer naming files explicitly instead of routinely using `git add .`. This
reduces the chance of committing generated results, large datasets, temporary
files, or credentials. If every displayed change is intentional, stage all of
them with `git add -A`, then check `git status` and `git diff --staged` before
committing.

## Bringing in collaboration updates

Fetch the latest state of the original repository:

```bash
git fetch upstream
```

Update the clean `master` branch using a fast-forward only:

```bash
git switch master
git merge --ff-only upstream/master
git push origin master
```

Return to the analysis branch and merge the updated infrastructure:

```bash
git switch analysis/hadronic-susy
git merge master
```

If Git reports conflicts, resolve each conflicted file, then stage the resolved
files and finish the merge:

```bash
git status
git add path/to/resolved-file
git commit
```

Run the test suite after merging:

```bash
.venv-local10k/bin/pytest -q
```

If the tests pass, push the updated analysis branch:

```bash
git push
```

This structure allows upstream infrastructure updates to enter the analysis
without requiring analysis files to be merged back into the original project.
Deleting or reorganizing files on `analysis/hadronic-susy` does not prevent
future pulls. Git records those changes, and only asks for manual conflict
resolution if a future upstream update edits the same parts of the same files.

## Useful inspection and recovery commands

Show the current branch and uncommitted changes:

```bash
git status
```

Show unstaged changes:

```bash
git diff
```

Show staged changes:

```bash
git diff --staged
```

Show recent commits and branch positions:

```bash
git log --oneline --decorate -10
```

Remove a file from the staging area without deleting its edits:

```bash
git restore --staged path/to/file
```

Restore uncommitted edits to a file only when those edits are definitely no
longer wanted:

```bash
git restore path/to/file
```

Do not use `git reset --hard` as a routine recovery command; it can permanently
discard work.

## Quick reference

For ordinary changes, the short version is:

```bash
git switch analysis/hadronic-susy
git status
git add path/to/changed-file
git diff --staged
git commit -m "Describe the change"
git push
```
