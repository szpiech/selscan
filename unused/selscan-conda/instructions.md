


```
conda build . --output-folder ./builds
````


anaconda upload builds/linux-64/selscan-1.0.0-0.tar.bz2
anaconda upload builds/osx-64/selscan-1.0.0-0.tar.bz2
anaconda upload builds/win-64/selscan-1.0.0-0.tar.bz2


conda install -c bioconda selscan



test
```
conda build . --no-test --output-folder build --platform win-64
```