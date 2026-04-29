#!/usr/bin/env sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
reference_dir=${LETS_BE_RATIONAL_CPP_DIR:-"$repo_root/../vollib_reference/LetsBeRational/LetsBeRational"}
out_dir="$repo_root/test/generated"
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/lbr-cpp-vectors.XXXXXX")
binary="$build_dir/generate_lets_be_rational_vectors"
json="$out_dir/LetsBeRationalCppVectors.json"
trap 'rm -rf "$build_dir"' EXIT

mkdir -p "$out_dir"

perl -Mutf8 -Mopen=:std,:encoding\(UTF-8\) -pe 's/ν/nu/g; s/h²/h_sqr/g; s/t²/t_sqr/g; s/h₂/hh2/g; s/h₃/hh3/g; s/h₄/hh4/g; s/q₁/q1/g; s/q₂/q2/g; s/θ/theta/g; s/𝛽/beta/g; s/sₗ/s_l/g; s/bₗ/b_l/g; s/μ/mu/g; s/λ/lambda/g; s/rₗₗ/r_ll/g; s/rₗₘ/r_lm/g; s/vₗ/v_l/g; s/sᵤ/s_u/g; s/bᵤ/b_u/g; s/rᵤₘ/r_um/g; s/vᵤ/v_u/g; s/rᵤᵤ/r_uu/g; s/x²/x_sqr/g; s/s³/s_cube/g;' \
  "$reference_dir/lets_be_rational.cpp" > "$build_dir/lets_be_rational.cpp"
perl -Mutf8 -Mopen=:std,:encoding\(UTF-8\) -pe 's/½/_half/g; s/𝑒/e/g; s/₀/0/g; s/ˣ/x/g; s/²/_sqr/g; s/√/sqrt/g; s/π/pi/g;' \
  "$reference_dir/normaldistribution.cpp" > "$build_dir/normaldistribution.cpp"

c++ -std=c++17 -O2 -DNO_XL_API \
  -I "$reference_dir" \
  "$repo_root/scripts/cpp/generate_lets_be_rational_vectors.cpp" \
  "$reference_dir/erf_cody.cpp" \
  "$build_dir/lets_be_rational.cpp" \
  "$build_dir/normaldistribution.cpp" \
  "$reference_dir/rationalcubic.cpp" \
  -o "$binary"

"$binary" > "$json"
printf '%s\n' "$json"
