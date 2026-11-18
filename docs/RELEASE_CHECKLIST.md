# Preparing or updating a release

The GitHub repository URL is recorded in `CITATION.cff`. An article DOI remains
omitted until one exists. The checklist below documents the initial publication
procedure and the steps to use for a later release.

1. Publish only this repository directory. Manuscripts, editorial letters and
   the author's submission instructions are distributed separately.
2. Commit with your own Git identity. Publish the repository under your chosen
   account, then check the Actions run. The workflow checks saved-output
   reconstruction; it does not run MATLAB or regenerate interval proofs.
3. Once the repository exists, add its actual URL as `repository-code` in
   `CITATION.cff`. Add a paper DOI only after one has been assigned.
4. After an intentional metadata edit, run, from the repository root:

   ```text
   python scripts/update_manifest.py
   python scripts/verify_release.py
   ```

   The manifest updater verifies the files listed in `SOURCE_MANIFEST.json`
   against their frozen hashes before rewriting `RELEASE_MANIFEST.json`.
   It excludes `.git/`, `generated/` and Python caches. It is not a substitute
   for validating a scientific change. Preserve the frozen archive when
   developing another study.
5. Commit the updated metadata and manifest, then create the GitHub release
   with tag `v1.0.0` and the description in `CHANGELOG.md`. If files have changed
   since the supplied ZIP was prepared, use a fresh archive of that commit;
   do not label the old ZIP as containing the new metadata.
6. Use the actual repository URL in the manuscript's Data Availability
   Statement. This local package does not claim that publication has occurred.

The journal's `Computational_Materials.zip` and the supplied GitHub ZIP have
the same file contents at handoff. Their ZIP directory prefixes differ. Keep
an archived copy of that handoff if later updating the public repository.
