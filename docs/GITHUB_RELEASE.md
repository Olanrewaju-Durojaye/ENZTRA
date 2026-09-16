# Publishing ENZTRA v1.0.1 on GitHub

## Repository contents

Commit the source tree, tests, documentation, examples, configuration example,
license, citation metadata, and GitHub Actions workflow. Do not commit local
environments, configured absolute paths, model weights, containers, or `jobs/`
directories; these are excluded by `.gitignore`.

## Release checklist

1. Install the archive in a clean environment.
2. Confirm `pytest` passes and `enztra.__version__` is `1.0.1`.
3. Confirm `examples/release_500/qualification_summary.json` matches the
   preserved full qualification report.
4. Create and push the `v1.0.1` annotated tag.
5. Create a GitHub release from that tag and attach the source ZIP.
6. Describe predictions as computational hypotheses requiring experimental
   validation and cite each independently installed upstream model.

Suggested commands after creating an empty GitHub repository:

```bash
git init
git add .
git commit -m "Release ENZTRA v1.0.1"
git branch -M main
git remote add origin <YOUR_GITHUB_REPOSITORY_URL>
git push -u origin main
git tag -a v1.0.1 -m "ENZTRA v1.0.1"
git push origin v1.0.1
```

Replace the placeholder with the repository URL shown by GitHub. Review staged
files with `git status` before committing; the generated 500-design job remains
local unless a deliberately curated, non-sensitive data deposit is prepared.
