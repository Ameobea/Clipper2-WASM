#!/usr/bin/env bash

# Set default build type to 'dev' if not specified
BUILD_TYPE=${1:-dev}

# Create necessary directories
mkdir -p clipper2-wasm/dist
mkdir -p clipper2-wasm/dist/es
mkdir -p clipper2-wasm/dist/umd

# Path to Clipper2 source (using local modified version, not submodule)
# Assumes this script is run from the Clipper2-WASM directory which is inside Clipper2/
CLIPPER2_PATH="../CPP"
CLIPPER2_SRC="${CLIPPER2_PATH}/Clipper2Lib/src"

# Clipper2 source files needed for offset functionality
# Note: We compile these directly with em++ rather than linking against a static lib
# because the static lib was compiled with GCC, not Emscripten
CLIPPER2_SOURCES="${CLIPPER2_SRC}/clipper.engine.cpp ${CLIPPER2_SRC}/clipper.offset.cpp ${CLIPPER2_SRC}/clipper.rectclip.cpp"

# Common flags
COMMON_FLAGS="-I${CLIPPER2_PATH}/Clipper2Lib/include -DUSINGZ --bind -s MODULARIZE=1 -s WASM_BIGINT -s ALLOW_MEMORY_GROWTH=1 -s EXIT_RUNTIME=0"

# Development build flags
DEV_FLAGS="-g3 -s ASSERTIONS=2 -s DISABLE_EXCEPTION_CATCHING=0 -s SAFE_HEAP=1 -O0"

# Production build flags
PROD_FLAGS="-O3 -s DISABLE_EXCEPTION_CATCHING=1"

# Select flags based on build type
if [ "$BUILD_TYPE" = "dev" ]; then
    FLAGS="$COMMON_FLAGS $DEV_FLAGS"
    echo "Building Development Version"
else
    FLAGS="$COMMON_FLAGS $PROD_FLAGS"
    echo "Building Production Version"
fi

# Check if Clipper2 source files exist
if [ ! -f "${CLIPPER2_SRC}/clipper.offset.cpp" ]; then
    echo "Error: Clipper2 source files not found at ${CLIPPER2_SRC}!"
    exit 1
fi

echo "Using Clipper2 source from: ${CLIPPER2_PATH}"

# Build ES6 module
echo "Building ES6 module..."
em++ $FLAGS \
    -s EXPORT_ES6=1 \
    -s NO_FILESYSTEM=1 \
    -s ENVIRONMENT='web' \
    -s EXPORT_NAME="Clipper2Z" \
    clipper2-wasm/clipper.bindings.cpp \
    ${CLIPPER2_SOURCES} \
    -o clipper2-wasm/dist/es/clipper2z.js \
    --post-js clipper2-wasm/glue-stub-z.js

if [ $? -eq 0 ]; then
    echo "ES6 build successful!"
else
    echo "ES6 build failed!"
    exit 1
fi

# Build UMD module
echo "Building UMD module..."
em++ $FLAGS \
    -s NO_FILESYSTEM=1 \
    -s ENVIRONMENT='web' \
    -s EXPORT_NAME="Clipper2ZFactory" \
    clipper2-wasm/clipper.bindings.cpp \
    ${CLIPPER2_SOURCES} \
    -o clipper2-wasm/dist/umd/clipper2z.js \
    --post-js clipper2-wasm/glue-stub-z.js

if [ $? -eq 0 ]; then
    echo "UMD build successful!"
else
    echo "UMD build failed!"
    exit 1
fi

echo ""
echo "Build complete!"
echo "Output files:"
ls -lh clipper2-wasm/dist/es/clipper2z.* 2>/dev/null
ls -lh clipper2-wasm/dist/umd/clipper2z.* 2>/dev/null

# ============================================================================
# COMMENTED OUT: Tools module build (SVG utilities - not needed for offset-only)
# ============================================================================
# CLIPPER2_UTILS_SRC="${CLIPPER2_PATH}/Utils"
# em++ $FLAGS -I${CLIPPER2_UTILS_SRC} -s EXPORT_ES6=1 -s ENVIRONMENT='web' clipper2-wasm/clipper-tools.bindings.cpp ${CLIPPER2_UTILS_SRC}/ClipFileSave.cpp ${CLIPPER2_UTILS_SRC}/ClipFileLoad.cpp ${CLIPPER2_UTILS_SRC}/clipper.svg.cpp ${CLIPPER2_SOURCES} -o clipper2-wasm/dist/es/clipper2z-utils.js -s EXPORT_NAME="Clipper2ZUtils" -s EXPORTED_FUNCTIONS=[FS] --post-js clipper2-wasm/glue-stub-tools-z.js
