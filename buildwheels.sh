#!/bin/bash

set -e
set -x

PYTHON_DIR=/opt/python
PYTHON_VERSION=3
PYTHON_BIN=( "$PYTHON_DIR/cp$PYTHON_VERSION"*/bin )
PYTHON_BIN_PROTO="${PYTHON_BIN[0]}"

# Update package manager
echo "=== Updating yum ==="
yum update -y

# Upgrade pip
echo "=== Upgrading pip ==="

for BIN in "${PYTHON_BIN[@]}"; do
  "$BIN/pip" install --upgrade pip
done

# Check if package needs to be published
if [[ -z "$BUILDWHEELS_SKIP_CHECK_VERSION" ]]; then
  echo "=== Checking Package Version ==="

  yum install -y bc

  PACKAGE_NAME="$("$PYTHON_BIN_PROTO/python" "$SRC_DIR/setup.py" --name)"
  PACKAGE_VERSION="$("$PYTHON_BIN_PROTO/python" "$SRC_DIR/setup.py" --version)"

  PYPI_PACKAGE_OUTDATED=false

  for BIN in "${PYTHON_BIN[@]}"; do
    # This tries to find the versions of the package available on PyPi in an
    # overly complicated manner (due to the fact that pip will not simply spit out
    # this information in a sane way)

    PYPI_PACKAGE_VERSIONS="$(
      "$BIN/pip" install "$PACKAGE_NAME==" 2>&1 1>/dev/null | \
       sed '1q;d' | sed 's/.*(from versions: \(.*\)).*/\1/'
    )"

    if [[ "$PYPI_PACKAGE_VERSIONS" = "none" ]]; then
      PYPI_PACKAGE_OUTDATED=true
      break
    fi;

    NEWEST_PYPI_PACKAGE_VERSION="$(echo "$PYPI_PACKAGE_VERSIONS" | sed 's/, /\n/g' | tail -n 1)"

    if [[ $(bc -l <<< "$PACKAGE_VERSION > $NEWEST_PYPI_PACKAGE_VERSION") -eq 1 ]]; then
      PYPI_PACKAGE_OUTDATED=true
      break
    fi
  done

  if [[ "$PYPI_PACKAGE_OUTDATED" = "false" ]]; then
    echo "PyPi package up-to-date, nothing to do"
    exit 0
  else
    echo "PyPi package outdated"
  fi
fi

# Install dependencies
if [[ -z "$BUILDWHEELS_SKIP_INSTALL_DEPS" ]]; then
  # Install Boost
  echo "=== Installing Boost ==="

  yum groupinstall -y "Development Tools"

  BOOST_DIR=/tmp/boost
  BOOST_VERSION=1.72.0
  BOOST_VERSION_="$(echo "$BOOST_VERSION" | tr . _)"
  BOOST_ARCHIVE="boost_$BOOST_VERSION_.tar.gz"

  if [[ ! -z "$BUILDWHEELS_DEBUG" ]]; then
    BOOST_VARIANT="variant=debug"
  else
    BOOST_VARIANT="variant=release"
  fi

  if [[ ! -z "$BUILDWHEELS_LINK_STATIC" ]]; then
    BOOST_LINK="link=static"
  else
    BOOST_LINK="link=shared"
  fi

  mkdir -p "$BOOST_DIR"
  cd "$BOOST_DIR"
  curl -L "https://sourceforge.net/projects/boost/files/boost/$BOOST_VERSION/$BOOST_ARCHIVE" -o "$BOOST_ARCHIVE"
  tar zxf "$BOOST_ARCHIVE"
  cd "boost_$BOOST_VERSION_"
  ./bootstrap.sh --with-libraries=graph
  ./b2 "$BOOST_VARIANT" "$BOOST_LINK"
  cd -

  # Install Lua
  echo "=== Installing Lua ==="

  yum install -y readline-devel

  LUA_DIR=/tmp/lua
  LUA_VERSION=5.2.0
  LUA_ARCHIVE="lua-$LUA_VERSION.tar.gz"

  if [[ ! -z "$BUILDWHEELS_DEBUG" ]]; then
    LUA_DEBUG_FLAGS="-Og -g"
  fi

  mkdir -p "$LUA_DIR"
  cd "$LUA_DIR"
  curl "https://www.lua.org/ftp/$LUA_ARCHIVE" -o "$LUA_ARCHIVE"
  tar zxf "$LUA_ARCHIVE"
  cd "lua-$LUA_VERSION"
  make linux MYCFLAGS="-fPIC $LUA_DEBUG_FLAGS"
  make install
  cd -

  # Install Cmake
  echo "=== Installing CMake ==="

  "$PYTHON_BIN_PROTO/pip" install cmake
  rm -f /usr/bin/cmake
  ln -s "$PYTHON_BIN_PROTO/cmake" /usr/bin/cmake

  # Install Twine
  echo "=== Installing Twine ==="
  "$PYTHON_BIN_PROTO/pip" install twine
fi

# Build wheels
if [[ -z "$BUILDWHEELS_SKIP_BUILD_WHEELS" ]]; then
  # Compile wheels
  echo "=== Compiling Wheels ==="

  export Boost_DIR="/tmp/boost/boost_1_72_0/stage/lib/cmake"

  if [[ ! -z "$BUILDWHEELS_DEBUG" ]]; then
    DEBUG_BUILD="debug-build=ON"
  fi

  for BIN in "${PYTHON_BIN[@]}"; do
    cd "$SRC_DIR"
    "${BIN}/python" "$SRC_DIR/setup.py" \
      build_ext "$DEBUG_BUILD" \
      bdist_wheel -d "$WHEELS_DIR"
    rm -rf *egg-info
    cd -
  done

  # Run auditwheel
  echo "=== Running auditwheel ==="

  export LD_LIBRARY_PATH="$BOOST_DIR/boost_$BOOST_VERSION_/stage/lib:$LD_LIBRARY_PATH"

  for WHEEL in "$WHEELS_DIR"/*.whl; do
    auditwheel repair "$WHEEL" -w "$WHEELS_DIR"
    rm "$WHEEL"
  done
fi

# Upload wheels
if [[ -z "$BUILDWHEELS_SKIP_UPLOAD_WHEELS" ]]; then
  echo "=== Uploading Wheels ==="

  for WHEEL in "$WHEELS_DIR"/*.whl; do
    "$PYTHON_BIN_PROTO/twine" upload \
        --skip-existing \
        -u "$TWINE_USER_NAME" \
        -p "$TWINE_PASSWORD" \
        "$WHEEL"
  done
fi
