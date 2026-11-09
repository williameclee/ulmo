# ULMO

## Cache handling

- Keep cache filenames stable. Never add a cache version to a filename or store or validate a cache version in the file.
- When reading a cache, check only for the variables the reader needs. Regenerate the data if any are missing and use the same path when saving.
- Leave clearing outdated caches to the user. Preserve each function's documented refresh and save options.
