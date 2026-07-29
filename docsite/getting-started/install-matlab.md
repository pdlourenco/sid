# Install — MATLAB / Octave

Pure MATLAB/Octave code with **zero toolbox dependencies**. Tested on
**MATLAB R2016b+** and **GNU Octave 8+** in CI on every commit.

## Clone and add to path

```bash
git clone https://github.com/pdlourenco/sid.git
```

In a MATLAB or Octave session:

```matlab
run('/path/to/sid/matlab/sidInstall.m')
```

To make the path persistent across sessions, add the line above to your
[`startup.m`](https://www.mathworks.com/help/matlab/ref/startup.html)
(or [`.octaverc`](https://docs.octave.org/latest/Startup-Files.html) for
Octave). No `pkg install` is needed — `sid` is a plain directory of `.m`
files that runs identically on both platforms.

## Quick smoke test

```matlab
help sidFreqBT
```

If the H1 line and parameter listing print, the path is wired up
correctly. Head to [Quick start](quick-start.md) for a runnable example,
or jump to the [MATLAB script catalog](../examples/matlab/index.md).

## Compatibility

| Platform | Minimum version |
|---|---|
| MATLAB | R2016b |
| GNU Octave | 8.0 |

`sid` deliberately avoids any Toolbox dependency, so the same scripts run
on a base MATLAB licence and on Octave.
