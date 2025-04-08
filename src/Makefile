# Detect OS
UNAME_S := $(shell uname -s)
ifeq ($(OS),Windows_NT)
    OS_TYPE := Windows
else ifeq ($(UNAME_S),Linux)
    OS_TYPE := Linux
else ifeq ($(UNAME_S),Darwin)
    OS_TYPE := Mac
else
    OS_TYPE := Unknown
endif

# Define default target
all:
ifeq ($(OS_TYPE),Linux)
	@$(MAKE) -f Makefile_linux
else ifeq ($(OS_TYPE),Mac)
	@$(MAKE) -f Makefile_macos
else ifeq ($(OS_TYPE),Windows)
	@$(MAKE) -f Makefile_win
else
	@echo "Unsupported OS"
endif

# Clean target for all platforms
clean:
ifeq ($(OS_TYPE),Linux)
	@$(MAKE) -f Makefile_linux clean
else ifeq ($(OS_TYPE),Mac)
	@$(MAKE) -f Makefile_macos clean
else ifeq ($(OS_TYPE),Windows)
	@$(MAKE) -f Makefile_win clean
else
	@echo "Unsupported OS: $(OS_TYPE)"
endif