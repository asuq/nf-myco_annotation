#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd -P)"
java_bin="java"
real_path="$(readlink -f "${script_dir}/minced.real")"
jar_dir="$(dirname "${real_path}")"

if [[ -n "${JAVA_HOME:-}" && -x "${JAVA_HOME}/bin/java" ]]; then
    java_bin="${JAVA_HOME}/bin/java"
fi

exec "${java_bin}" -XX:+PerfDisableSharedMem -jar "${jar_dir}/minced.jar" "$@"
