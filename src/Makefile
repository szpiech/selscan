# Makefile for selscan: Linux, macOS, and Windows (MSYS2 / MinGW-w64).
#
#   make            build ./selscan
#   make check      build, then run the regression suite in ../tests
#   make install    copy selscan to $(PREFIX)/bin (PREFIX defaults to /usr/local)
#   make clean      remove the build directory and the binary
#
# Needs a C++17 compiler and zlib (headers and library). Set any of these
# on the command line or in the environment:
#
#   CXX             compiler (make's default is g++, which is clang on macOS)
#   CXXFLAGS        optimisation and debugging flags (default -O3)
#   CPPFLAGS, LDFLAGS, LDLIBS
#   STATIC=1        Linux: link libstdc++, libgcc and zlib statically and
#                   keep glibc dynamic, as the release binaries are built.
#                   Windows: a fully static .exe.
#   UNIVERSAL=1     macOS: one binary for arm64 and x86_64.
#   MACOS_MIN=11.0  macOS: oldest macOS the binary runs on.
#
# Example debug build with UBSan:
#   make CXXFLAGS="-O1 -g -fsanitize=undefined" LDFLAGS=-fsanitize=undefined

CXXFLAGS ?= -O3
override CXXFLAGS += -std=c++17 -Wall -Wextra
override CPPFLAGS += -MMD -MP
PREFIX ?= /usr/local

ifeq ($(OS),Windows_NT)
  PLATFORM := windows
else
  PLATFORM := $(shell uname -s | tr '[:upper:]' '[:lower:]')
endif

EXE :=
SYSLIBS := -lz
ifeq ($(PLATFORM),windows)
  EXE := .exe
  ifeq ($(STATIC),1)
    override LDFLAGS += -static
  endif
else ifeq ($(PLATFORM),darwin)
  ifeq ($(UNIVERSAL),1)
    override CXXFLAGS += -arch arm64 -arch x86_64
  endif
  ifneq ($(MACOS_MIN),)
    override CXXFLAGS += -mmacosx-version-min=$(MACOS_MIN)
  endif
else
  ifeq ($(STATIC),1)
    override LDFLAGS += -static-libstdc++ -static-libgcc
    SYSLIBS := -Wl,-Bstatic -lz -Wl,-Bdynamic
  endif
  SYSLIBS += -lm
endif
override LDLIBS += $(SYSLIBS) -pthread

SRCS := selscan-main.cpp selscan-cli.cpp selscan-data.cpp param_t.cpp \
        norm.cpp gene.cpp binom.cpp gzstream.cpp thread_pool.cpp \
        hapmap/hapdata.cpp hapmap/mapdata.cpp hapmap/bitset.cpp \
        stats/ehh.cpp stats/ehh12.cpp stats/ihs.cpp stats/ihh12.cpp \
        stats/pi.cpp stats/xpihh.cpp
BUILD := build
OBJS := $(SRCS:%.cpp=$(BUILD)/%.o)
BIN := selscan$(EXE)

all: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# header dependencies written by -MMD
-include $(OBJS:.o=.d)

check: $(BIN)
	../tests/run_tests.sh --bin ./$(BIN)

install: $(BIN)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/

clean:
	rm -rf $(BUILD) $(BIN)

.PHONY: all check install clean
