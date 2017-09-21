Small library of open source VTK code that was started in the TU Delft
visualization group.

This library was resurrected and partially ported to VTK 7.1 in 2017.

# Installing from anaconda

```
conda config --add channels cpbotha
conda install vtktudoss
```

Alternatively:

```
conda install -c cpbotha vtktudoss
```

# Building conda packages from source

```
cd vtktudoss/conda.recipe
conda build .
```

# Notes

The wrapping setup in `Contrib/STLib` is the most up to date at the
moment.
