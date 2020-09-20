# LSynth

LDraw CAD system compatible flexible part synthesizer.  This program
reads your LDraw file, searching for unofficial META commands (structured
comments), that specify things you want synthesized.

### What does it do?

LSynth has two primary forms of synthesis.  The first kind of synthesis
creates hose like things including:
    rubbed, rubber, pneumatic, and flex system hoses, as well as electric
    and fiber optic cables, flex system cables, flexible axles, string, and
    minifig chain.

It also creates things that travel around circular lego parts.  Things
like:
    rubber band, rubber belt, technic chain, technic plastic tread, technic
    rubber tread.

The files `tube.c`, `tube.h`, `curve.c` and `curve.h` perform hose synthesis.
The files `band.c` and `band.h` perform band synthesis.

This `lsynthcp.c` file contains the main entry/exit points for the program.
It opens and scans the LDraw file provided, identifies synthesis
synthesis specifications and hands them off to the appropriate synthesis
methodology.

### Obtaining LSynth binaries

You may be able to obtain pre-built binaries from [here](http://lsynth.sourceforge.net/).
Place the executable in this directory.

To build it from source instead, just run:

```
make
```

### Usage

`lsynthcp [-v] [-h] [-m] [-l] [-p] <src> <dst>`

- -v: prints `lsynthcp` version
- -h: prints this help message
- -m: prints out the LSynth portion of the `MLcad.ini` for using this program
- -l: format the output as an official ldraw part
- -p: prints out the full path name of the this executable

See [Willy Tscager's tutorial
page](http://www.holly-wood.it/lsynth-en.html) for using LSynth within MLCAD.

To create a flexible part, you put specifications for the part
directly into your LDraw file, where the part is needed.

Example usage to recreate one of the provided example files:

```
./lsynthcp examples/ELECTRIC_NXT_CABLE-Constraints.ldr examples/result.ldr
```
