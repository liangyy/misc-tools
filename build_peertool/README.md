As I encountered many failures on building `peertool`, I left this note to record the occasional success for future reference.

# Prerequisite

Compiler: `gcc-5`
Python: `python2.7`

# Download the source

```
git clone https://github.com/PMBio/peer
```

Or download the zip file.

# Build

```
cd peer
mkdir build; cd build
cmake -D BUILD_PEERTOOL=1  -DCMAKE_INSTALL_PREFIX=path-to-install ./..
make
make install
```

# Result

The executable `peertool` will be installed at `path-to-install/bin`
