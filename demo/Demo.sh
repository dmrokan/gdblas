#!/bin/sh
echo -ne '\033c\033]0;Demo\a'
base_path="$(dirname "$(realpath "$0")")"
"$base_path/Demo.x86_64" "$@"
