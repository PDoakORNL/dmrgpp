# dmrgpp / cincuenta — Claude Code notes

## Temp directory

`/private/tmp` fills up frequently on this machine. Always use the project-local tmp
directory for any shell commands that need a writable scratch area:

```
/Users/Shared/ornldev/code/dmrgpp/build/tmp
```

Set this as `CLAUDE_CODE_TMPDIR` when needed, or redirect output there directly.
Never rely on `/tmp` or `/private/tmp`.

## Build layout

- Source root: `/Users/Shared/ornldev/code/dmrgpp`
- Build tree:  `/Users/Shared/ornldev/code/dmrgpp/build`
- cincuenta executable: `build/cincuenta/src/cincuenta`
- Run tests from: `build/cincuenta/src/` with `ctest`

## Running tests

```bash
cd /Users/Shared/ornldev/code/dmrgpp/build/cincuenta/src
ctest -LE Nightly          # fast CI suite (excludes slow Nightly tests)
ctest -R <pattern> -V      # run specific test with verbose output
```

Build first: `cmake --build /Users/Shared/ornldev/code/dmrgpp/build --target cincuenta -j4`

Note: `build/cincuenta/src/CTestTestfile.cmake` is cmake-generated. When changing test
reference values, update **both** `cincuenta/src/CMakeLists.txt` (authoritative source)
and `build/cincuenta/src/CTestTestfile.cmake` (used by ctest until next cmake regeneration).
