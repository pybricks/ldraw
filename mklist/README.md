# mklist

mklist v1.6 20100706 (C) 1999-2010 Lars C. Hassing. Replacement for James Jessiman's makelist

### What does it do?

It parses the `ldraw/parts` folder to produce a text file with all available parts.

### Obtaining mklist binaries

The included [`mklist.exe`](../mklist.exe) file is obtained from the official LDraw
archives, for full compatibility with the official parts archive.

To build it from source instead, run:

```
make
```

### Usage

`mklist options`

Options:
- -h: You already figured this one out :-)
- -n: Sort by Number
- -d: Sort by Description
- -c: Check for duplicate descriptions. "parts.lst" unchanged.
- -m: Don't skip parts with "~Moved to xxx" description
- -~: Skip parts with ~ description, e.g. "~Winch  2 x  4 x  2 Top"
- -i `<dir>`: input directory, default is "PARTS" in current directory
- -o `<file>`: output filename, default is "parts.lst" in current directory
- -f: Force it to finish.  No prompts.
- -q: Quiet mode.  No warnings, and no prompts.
- -8: Use 8.3 names for compatibility.
- -t: Truncating descriptions to fit in an 80 char terminal window.
- -r: Ragged filename column.  Size it to fit short filenames.
- -l: Truncate Long descriptions at 64 chars.
- -v: Print verbose info.  Useful for debugging.

Typical usage to produce the usual `parts.lst` in the `ldraw` root folder:

```
./mklist -d -i ../parts -o ../parts.lst
```
