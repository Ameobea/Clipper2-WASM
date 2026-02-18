# Clipper2-WASM build tasks

# Build the WASM module (production by default)
build mode="prod":
  #!/bin/bash
  source ~/emsdk/emsdk_env.sh 2>/dev/null
  bash clipper2-wasm/compile-wasm.sh {{mode}}

  # Optimize WASM binary
  if command -v wasm-opt &> /dev/null; then
    echo "Optimizing WASM with wasm-opt..."
    wasm-opt clipper2-wasm/dist/es/clipper2z.wasm -g -c -O4 \
      --enable-bulk-memory \
      --enable-nontrapping-float-to-int \
      --precompute-propagate \
      --strip-dwarf \
      -o clipper2-wasm/dist/es/clipper2z.wasm
    echo "Optimization complete!"
    ls -lh clipper2-wasm/dist/es/clipper2z.wasm
  else
    echo "wasm-opt not found, skipping optimization"
  fi

# Build development version with debug symbols
build-dev:
  just build dev

# Install ES6 module to dream project
install:
  mkdir -p ~/dream/src/viz/wasm/clipper2
  cp clipper2-wasm/dist/es/clipper2z.js ~/dream/src/viz/wasm/clipper2/
  cp clipper2-wasm/dist/es/clipper2z.wasm ~/dream/src/viz/wasm/clipper2/
  echo "Installed to ~/dream/src/viz/wasm/clipper2/"
  ls -lh ~/dream/src/viz/wasm/clipper2/

# Build and install in one step
deploy: build install

# Clean build artifacts
clean:
  rm -rf clipper2-wasm/dist/es/*
  rm -rf clipper2-wasm/dist/umd/*
  echo "Cleaned build artifacts"
