Basis: [Clean Code](https://gist.github.com/wojteklu/73c6914cc446146b8b533c0988cf8d29) and [PEP 8](https://peps.python.org/pep-0008/) but not up to the comma.

- General:
    - Code is simple and explains itself without comments or documentation.
    - Large constructs are broken down into small ones with a single responsibility. Each one should fit one computer screen.
    - Abstraction layers are well separated.

- Procedures: subroutines and functions
    - `pure` without side effects whenever possible.
    - Avoid long argument lists by introducing structs and/or classes.
    - bundle together into a static `module` (no `module` variables)

- Structs
    - To combine the same data that are used by many procedures.
    - `public` data, procedures live outside
    - `type` as struct: preferred
        - Multiple instances can exist.

- Classes
    - To decouple abstract interface from concrete implementations and/or ...
    - ... to combine many procedures that use the same data.
    - `public` subroutines and functions, `protected` (syntax or convention) data.
    - `type` as class: preferred
        - multiple instances can exist
        - derive from abstract type and use factory pattern

- Arrays
    - are to be indexed starting by 1

- Module as singleton ...
    - ... struct: legacy
        - as single instance in global namespace
        - avoid state changes wherever possible
        - put procedures using the variables in paired `modulename_sub` for safety (see `module` naming)
    - ... class: in justified cases
        - as single instance in global namespace
        - only do small, single purpose `modules`
        - can contain procedures and mutable `module` variables
        - e.g. config `module`: load config from input file and save in `read-only` variables

- Tests:
    - Write test *before* code
    - In classes, prefer testing of public methods to state-based testing
    - tests are located in folders `test`

- Documentation:
    - Prefer code that speaks for itself and scripts to docs
    - Save reusable code or script *before* typing in interactive shell
    - Priority for user docs: Examples > Quickstart > API docs
    - API docs may be made with the doxygen syntax (`!>` and `!<`) but in *plain text*
    - API docs comments may be made before function, type, module, ... declaration (headline-like)

- Git:
    - Commit message of maximum 50 characters in first line. Two styles possible
        - General comment starting uppercase and without final dot e.g.
            - `Change core logics in vmec_to_boozer.py`
        - Special prefixes, all lowercase: (`wip, clean, refactor, doc`) e.g.
            - `clean: remove unused code`
        - If required, add blank line and details below
    - Pull request title should specify which issue it closes, e.g.
        - `Fix #1234: Fortran GUI too beautiful`
    - Stale branches should be deleted or if they contain useful information, renamed to `archive/...` with an appropriate README

- Style:
    - Indent with four spaces, no tabs.
    - Max line width: 88 (like Python black)
    - `modules`
        - Fortran:
            - indent after module declaration
            - do `use modulename, only:`, do not `use modulename` for import
        - Python:
            - no wildcards for import (only if import few and small `modules`)
    - procedures
        - indent after procedure header
        - Prefer verbs as names, with exceptions such as `from_file()` returning a value.
        - indent as level marker -> all procedures with same indent (black)
        - 1 blank line between procedures

- Numerics
    - specify `real64` for floats
      ```fortran
      use iso_fortran_env, only: dp => real64
      real(dp) :: float_number
      ```

- Always declare `implicit none` (Fortran) ...
    - ... once in each `module` declaration
    - ... once in the `program` declaration

- Directories and file names:
    - lowercase with underlines
    - exceptions are for external standards like `CMake`

- `module` names:
    - preferred: `project_modulename` e.g. `libneo_utils`
      ```fortran
      module libneo_utils
          contains
          subroutine find_local_maximum(func)
          subroutine find_local_minimum(func)
      ```
        - special: `module` for `abstract type` (add key-suffix `_base`)
            ```fortran
            module libneo_field_base
                abstract type field_t
            ```
    - legacy ("`modules` as structs"):
        - `module` `modulename_mod` containing (mostly) only variables e.g. `arnoldi_mod`
        - `module` `modulename_sub` containing the procedures using these variables e.g. `arnoldi_sub`

- File i/o
    - read with `open(newunit=iunit, file=trim(filename), status='old', action='read')` (otherwise empty files are generated if they are missing)
    - write with `open(newunit=iunit, file=trim(filename), status='replace', action='write')`

