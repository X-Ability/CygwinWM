# Building MACE-enabled LAMMPS on Windows (MSYS2)

## 0. Prerequisites and Policy (Important)

- This targets the **CPU version only**. Because the MSYS2 LibTorch is a CPU build, the
  native Windows CUDA/GPU version is not covered here.
- **Use the UCRT64 environment, not the existing MINGW64 environment.**
  The MSYS2 package `python-pytorch`, which includes LibTorch, exists only in the UCRT64
  (and CLANG64) environments and is not available in the MINGW64 environment used by the
  current README.

### Conditions regarding Python (assumptions of this procedure)

- **The generated `lmp.exe` does not link against Python at all** (no Python needed at runtime).
- MSYS2 has no standalone LibTorch package, so you need to install
  `mingw-w64-ucrt-x86_64-python-pytorch`, which includes the C++ headers, DLLs, and `TorchConfig.cmake`.
- **Model conversion** (`create_lammps_model.py`) requires Python, but this is done in advance on a
  separate machine/environment to prepare the TorchScript model.

---

## 1. Installing the UCRT64 Toolchain and Required Packages

Launch the **"MSYS2 UCRT64"** shortcut (`ucrt64.exe`) and run the following.
(Note: this is *not* the MINGW64 shortcut.)

```shell
# Update the package database
$ pacman -Sy
$ pacman -Su

# Build tools
$ pacman -S --needed \
    mingw-w64-ucrt-x86_64-gcc \
    mingw-w64-ucrt-x86_64-gcc-fortran \
    mingw-w64-ucrt-x86_64-cmake \
    mingw-w64-ucrt-x86_64-ninja \
    mingw-w64-ucrt-x86_64-msmpi \
    git patch m4 make

# LibTorch
$ pacman -S --needed mingw-w64-ucrt-x86_64-python-pytorch
```

---
## 2. Prepare the LibTorch Import Libraries (Required)

The MSYS2 LibTorch ships **only the DLLs** for `libc10.dll` / `libtorch.dll` / `libtorch_cpu.dll`
and does not include the corresponding import libraries (`.dll.a`). However, the bundled
`Caffe2Targets.cmake` references these `.dll.a` files, so configuration stops with the following error:

```
CMake Error at .../share/cmake/Caffe2/Caffe2Targets.cmake:124 (message):
  The imported target "c10" references the file
     ".../torch/lib/libc10.dll.a"
  but this file does not exist.
```

Because GNU `ld` can take a PE-format DLL directly as a link input, **simply copying the three DLLs
under the `.dll.a` names** satisfies both the existence check and the linking:

```shell
$ cd /c/msys64/ucrt64/lib/python*/site-packages/torch/lib
$ cp libc10.dll libc10.dll.a
$ cp libtorch.dll libtorch.dll.a
$ cp libtorch_cpu.dll libtorch_cpu.dll.a
```

---

## 3. Obtaining the Source Code

```shell
$ cd ~
$ git clone --branch=mace --depth=1 https://github.com/ACEsuit/lammps.git lammps-mace
$ cd lammps-mace
```

---

## 4. Building with CMake

```shell
$ mkdir build && cd build
$ TORCH_ROOT=$(cygpath -m "$(ls -d /ucrt64/lib/python*/site-packages/torch)")
$ cmake -G Ninja ../cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_INSTALL_PREFIX="$PWD/../dist" \
    -D CMAKE_CXX_STANDARD=17 \
    -D CMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -D CMAKE_PREFIX_PATH="$TORCH_ROOT" \
    -D BUILD_MPI=ON \
    -D BUILD_OMP=ON \
    -D PKG_ML-MACE=ON \
    -D PKG_EXTRA-PAIR=ON \
    -D PKG_MOLECULE=ON \
    -D PKG_MISC=ON \
    -D PKG_RIGID=ON \
    -D PKG_REAXFF=ON \
    -D PKG_EXTRA-DUMP=ON \
    -D PKG_MC=ON \
    -D PKG_KSPACE=ON \
    -D PKG_PHONON=ON \
    -D PKG_QEQ=ON \
    -D PKG_REACTION=ON \
    -D PKG_REPLICA=ON \
    -D PKG_DPD-BASIC=ON \
    -D PKG_EXTRA-MOLECULE=ON \
    -D PKG_MANYBODY=ON \
    -D PKG_MEAM=ON \
    -D PKG_EXTRA-FIX=ON \
    -D PKG_FEP=ON \
    -D PKG_TALLY=ON \
    -D PKG_COLVARS=ON

$ cmake --build . -j
$ cmake --install .
```

On success, `/dist/bin/lmp.exe` is generated.

---

## 5. Packaging

`lmp.exe` depends on LibTorch and the MinGW runtime DLLs. Gather the DLLs required to run it.
**Do this in the UCRT64 shell** (so that `ldd` resolves the UCRT64 DLLs correctly).

```shell
$ cd ~/lammps-mace/dist/bin
$ strip lmp.exe         # strip the executable only (do not strip the DLLs)

# (a) Collect the full set of LibTorch runtime DLLs
#     (copy them all so that dynamically loaded ones are not missed)
$ cp -v /c/msys64/ucrt64/lib/python*/site-packages/torch/lib/*.dll .

# (b) Collect the UCRT64 MinGW runtime DLLs that the exe depends on
$ for dll in $(ldd lmp.exe | grep -i '/ucrt64/' | awk '{print $3}'); do cp -v "$dll" .; done
```

# (c) Add the MS-MPI runtime DLL
# Copy msmpi.dll from the LAMMPS distribution bundled with CygwinWM
# into this directory so that it is placed alongside lmp.exe.

Make sure that msmpi.dll is included in the final package together with lmp.exe.