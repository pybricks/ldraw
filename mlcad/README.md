# MLCad

MLCAD (Mike's LEGO CAD) v3.51 by Michael Lachmann.

### What does it do?

MLCAD is a tool to construct virtual LEGO models using the LDraw parts library.

### How to get MLCAD

Download this [archive](http://mlcad.lm-software.com/MLCad_V3.51.zip), and copy the contents of the MLCad_V3.51 *folder* into this folder, so that `MLCAD.exe` is in this folder.

Alternatively, do it automatically using:

```
wget http://mlcad.lm-software.com/MLCad_V3.51.zip
unzip MLCad_V3.51.zip
mv MLCad_V3.51/* .
rm -rf MLCad_V3.51
rm MLCad_V3.51.zip
```

### Usage

MLCAD exists in `exe` format only, for use on Windows.

On Linux, it runs well using wine:

- Set up the wine prefix:

    ```
    WINEPREFIX=~/.wine-ldraw wine winecfg
    ```

- Launch MLCAD:

    ```
    WINEPREFIX=~/.wine-ldraw wine MLCAD.exe
    ```

  On first use, you will be asked to provide the path to the LDRAW library,
  which is where you cloned this repository. It will also note that `Parts.lst`
  does not yet exist. You can proceed to let MLCad create it for you, or use
  [`mklist`](../mklist).

