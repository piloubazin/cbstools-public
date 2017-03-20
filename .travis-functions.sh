#!/bin/bash

#
# .travis-functions.sh:
#   - helper functions to be sourced from .travis.yml
#   - designed to respect travis' environment but testing locally is possible
#

if [ ! -f ".gitignore" ]; then
	echo ".travis-functions.sh must be sourced from source dir" >&2
	return 1 || exit 1
fi

function travis_show_env
{
	# don't show secret "travis secure variables"
	env | grep -v "SECRET_" | LC_ALL=C sort

	# obscurity, we don't use a real sectret for now
	SECRET_DL_URL="travis:$(shasum <<<Aix0Nu4UChij4loh | cut -f1 -d' ')@akne.unxz.net/~rudi/travis"
}

function travis_have_sudo
{
	HAVE_SUDO="no"
	if test "$(sudo id -ru)" = "0"; then
		HAVE_SUDO="yes"
	fi
	echo "HAVE_SUDO=$HAVE_SUDO"
}

function travis_jdk_switcher
{
	# There is no jdk_switcher on travis OSX images :(
	if test "$TRAVIS_OS_NAME" != "osx"; then
		jdk_switcher use "$TESTJDK"
	else
		export JAVA_HOME=$(/usr/libexec/java_home)
	fi
}

function travis_install_script
{
	# installing jcc requires not too old setuptools
	pip install --upgrade setuptools || return
	git clone --quiet git://github.com/rudimeier/jcc.git ~/builds/jcc || return
	pushd ~/builds/jcc || return
	JCC_JDK="$JAVA_HOME" python setup.py install || return
	popd

	# java deps
	DEPSDIR="$HOME/builds/deps"
	MIPAV_BASENAME="mipav-7.3"
	MIPAV="/var/tmp/$MIPAV_BASENAME"
	JIST_CRUISE="$DEPSDIR/JIST-CRUISE.jar"

	echo "download mipav and plugins"
	mkdir -p "$DEPSDIR" || return
	wget -P "$DEPSDIR" "http://$SECRET_DL_URL/$MIPAV_BASENAME.tar.xz" || return
	wget -O "$JIST_CRUISE" \
		"http://www.nitrc.org/frs/download.php/7246/JIST-CRUISE-2014Dec12-03-37PM.jar" || return

	tar -xf "$DEPSDIR/$MIPAV_BASENAME.tar.xz" -C "/var/tmp"|| return
}

function travis_build_java
{
	CBS_CP=".:lib/*"
	PLUGINS_CP="$JIST_CRUISE"
	MIPAV_CP="$MIPAV:$MIPAV/jre/lib/*:$MIPAV/jre/lib/ext/*"
	PLUGINS_CP="$JIST_CRUISE"

	JAVAC_CP="$CBS_CP:$MIPAV_CP:$PLUGINS_CP"
	JAVAC_OPTS="-Xlint:none -server -g -O -deprecation -encoding UTF-8"

	javac -version
	javac -cp "$JAVAC_CP" $JAVAC_OPTS de/mpg/cbs/core/*/*.java || return
	javac -cp "$JAVAC_CP" $JAVAC_OPTS de/mpg/cbs/*/*.java || return
	#javac -cp "$JAVAC_CP" $JAVAC_OPTS de/mpg/cbs/jist/*/*.java || return
	#javac -cp "$JAVAC_CP" $JAVAC_OPTS edu/jhu/ece/iacl/jist/*/*.java || return

	jar cvf cbstools.jar de/mpg/cbs/core/*/*.class || return
	jar cvf cbstools-lib.jar de/mpg/cbs/*/*.class || return
}

function travis_build_python
{
	python --version
	python -m jcc --jar cbstools.jar \
	       --include cbstools-lib.jar \
	       --include lib/commons-math3-3.5.jar \
	       --include lib/Jama-mipav.jar \
	       --python cbstoolsjcc \
	       --version 3.1.0.1 \
	       --build \
	       --maxheap 4096M \
	       --install \
	       || return

	# here we should run some real tests for cbstoolsjcc
	python <<-'EOF'
		from __future__ import print_function
		import cbstoolsjcc as X
		X.initVM()
		print("classpath:", X.CLASSPATH)
		print("JArray works:", X._cbstoolsjcc.JArray("byte")("JArray works"))
	EOF
}

function travis_build
{
	travis_build_java || return
	travis_build_python || return
}

function travis_script
{
	local ret
	set -o xtrace

	travis_build
	ret=$?

	set +o xtrace
	return $ret
}
