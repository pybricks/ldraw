# LDraw

This is a mirror of the LDraw LEGO parts library available at [ldraw.org](https://www.ldraw.org/). It is the same collection of parts, but with all the benefits of version control. See also the original LDraw [Readme.txt](Readme.txt).

### Usage

```
git clone https://github.com/pybricks/ldraw.git
```

This creates the `ldraw` folder. Use it like you always do.


### Changes compared to the official `complete.zip` parts pack:
- Full Git history since the 2008 re-release.
- Fixed [`mklist`](./mklist) source code. The `exe` files are unchanged.
- Includes [`lsynth`](./lsynth) source code and dependencies.

### Tips for Linux users (Optional)

The main purpose of this repository is simply to provide the LDraw library,
which is platform independent. Just clone the repository, or download and unzip [it](https://github.com/pybricks/ldraw/archive/master.zip).
You can do this on Windows, Mac, or Linux.

Most LDraw tools and tutorials assume that you use Windows. The following tips are intended to make life a bit easier for Linux users:
- Do not install `ldraw-parts` with `apt`. It is outdated.
- LSynth does not provide binaries for Linux, but you can [build it easily](./lsynth). All its dependencies are included, too. You no longer need to search the web for three different zip archives in different locations.
- A few Linux-compatible LDraw CAD tools do exist, such as [LeoCAD](https://www.leocad.org/). Also, [MLCAD works well with Wine.](./mlcad)
- [LPub3D](https://trevorsandy.github.io/lpub3d/) is available as an `AppImage`, which conveniently bundles several renderers with it. Before you run it, create a symlink to this repository using `ln -s my/path/to/repo/ldraw/ ~/Documents/ldraw` (replace the first path with your own). LPub3D will look there when you run it for the first time.

