#!/bin/bash
mkdir -p $PREFIX/bin

# Platform-specific copy
if [[ "$target_platform" == "linux-64" ]]; then
    cp bin/linux/selscan $PREFIX/bin/selscan
    cp bin/linux/norm $PREFIX/bin/norm
elif [[ "$target_platform" == "osx-64" ]]; then
    cp bin/osx/selscan $PREFIX/bin/selscan
    cp bin/osx/norm $PREFIX/bin/norm
elif [[ "$target_platform" == "osx-arm64" ]]; then
    cp bin/macos-arm64/selscan $PREFIX/bin/selscan
    cp bin/macos-arm64/norm $PREFIX/bin/norm
elif [[ "$target_platform" == "linux-aarch64" ]]; then
    cp bin/linux/selscan $PREFIX/bin/selscan
    cp bin/linux/norm $PREFIX/bin/norm
else
    echo "Unsupported platform: $target_platform"
    exit 1
fi

chmod +x $PREFIX/bin/selscan
