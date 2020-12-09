#!/bin/bash

set -e
set -x

PYTHON_DIR=/opt/python

# Locate Python versions >= 3.6
for BIN in "$PYTHON_DIR"/cp3*/bin; do
   MINOR_VERSION=$(echo "$BIN" | sed 's/.*cp3\([0-9]\+\).*/\1/')

   (( $MINOR_VERSION >= 6 )) && PYTHON_BIN+=("$BIN")
done

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
if [[ "$BUILDWHEELS_SKIP_CHECK_VERSION" != "y" ]]; then
  echo "=== Checking Package Version ==="

  yum install -y jq bc

  PACKAGE_NAME="$("$PYTHON_BIN_PROTO/python" "$SRC_DIR/setup.py" --name)"
  PACKAGE_VERSION="$("$PYTHON_BIN_PROTO/python" "$SRC_DIR/setup.py" --version)"

  PACKAGE_OUTDATED=false

  for BIN in "${PYTHON_BIN[@]}"; do
    # This tries to find the versions of the package available on PyPi in an
    # overly complicated manner (due to the fact that pip will not simply spit out
    # this information in a sane way)

    PYPI_JSON_URL="https://pypi.org/pypi/$PACKAGE_NAME/json"

    NEWEST_PACKAGE_VERSION=$(curl -L -s "$PYPI_JSON_URL" | \
                             jq  -r '.releases | keys | .[]' | \
                             sort -V | tail -n 1)

    if [[ $(bc -l <<< "$PACKAGE_VERSION > $NEWEST_PACKAGE_VERSION") -eq 1 ]]; then
      PACKAGE_OUTDATED=true
      break
    fi
  done

  if [[ "$PACKAGE_OUTDATED" = "false" ]]; then
    echo "PyPi package up-to-date, nothing to do"
    exit 0
  else
    echo "PyPi package outdated"
  fi
fi

# Install dependencies
if [[ "$BUILDWHEELS_SKIP_INSTALL_DEPS" != "y" ]]; then
  # Install Boost
  echo "=== Installing Boost ==="

  yum groupinstall -y "Development Tools"

  BOOST_DIR=/tmp/boost
  BOOST_VERSION=1.72.0
  BOOST_VERSION_="$(echo "$BOOST_VERSION" | tr . _)"
  BOOST_ARCHIVE="boost_$BOOST_VERSION_.tar.gz"

  if [[ ! -z "$BUILDWHEELS_DEBUG" ]]; then
    BOOST_VARIANT="variant=debug"
    BOOST_CXXFLAGS="-O2 -g3 -fno-omit-frame-pointer"
  else
    BOOST_VARIANT="variant=release"
    BOOST_CXXFLAGS="-O2"
  fi

  if [[ ! -z "$BUILDWHEELS_LINK_DYNAMIC" ]]; then
    BOOST_LINK="link=shared"
  else
    BOOST_LINK="link=static"
  fi

  mkdir -p "$BOOST_DIR"
  cd "$BOOST_DIR"
  curl -L "https://sourceforge.net/projects/boost/files/boost/$BOOST_VERSION/$BOOST_ARCHIVE" -o "$BOOST_ARCHIVE"
  tar zxf "$BOOST_ARCHIVE"
  cd "boost_$BOOST_VERSION_"
  ./bootstrap.sh --with-libraries=graph
  ./b2 "$BOOST_VARIANT" "$BOOST_LINK" cxxflags="$BOOST_CXXFLAGS"
  cd -

  # Install Lua
  echo "=== Installing Lua ==="

  yum install -y readline-devel

  LUA_DIR=/tmp/lua
  LUA_VERSION=5.2.0
  LUA_ARCHIVE="lua-$LUA_VERSION.tar.gz"

  if [[ ! -z "$BUILDWHEELS_DEBUG" ]]; then
    LUA_MYCFLAGS="-fPIC -O2 -g3 -fno-omit-frame-pointer"
  else
    LUA_MYCFLAGS="-fPIC -O2"
  fi

  mkdir -p "$LUA_DIR"
  cd "$LUA_DIR"
  curl "https://www.lua.org/ftp/$LUA_ARCHIVE" -o "$LUA_ARCHIVE"
  tar zxf "$LUA_ARCHIVE"
  cd "lua-$LUA_VERSION"
  make linux MYCFLAGS="$LUA_MYCFLAGS"
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
if [[ "$BUILDWHEELS_SKIP_BUILD_WHEELS" != "y" ]]; then
  # Compile wheels
  echo "=== Compiling Wheels ==="

  export Boost_DIR="/tmp/boost/boost_$BOOST_VERSION_/stage/lib/cmake"

  if [[ ! -z "$BUILDWHEELS_DEBUG" ]]; then
    BUILD_EXT="build_ext debug-build=ON"
  else
    BUILD_EXT="build_ext"
  fi

  for BIN in "${PYTHON_BIN[@]}"; do
    cd "$SRC_DIR"
    "${BIN}/python" "$SRC_DIR/setup.py" \
      "$BUILD_EXT" \
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
if [[ "$BUILDWHEELS_SKIP_UPLOAD_WHEELS" != "y" ]]; then
  echo "=== Uploading Wheels ==="

  for WHEEL in "$WHEELS_DIR"/*.whl; do
    "$PYTHON_BIN_PROTO/twine" upload \
        --repository "$TWINE_REPOSITORY" \
        --skip-existing \
        -u "$TWINE_USER_NAME" \
        -p "$TWINE_PASSWORD" \
        "$WHEEL"
  done
fi
