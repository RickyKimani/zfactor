# zfactor — common tasks.
#
# Run `just` with no arguments to see the list.

# Directory holding the page ranges the parsers were run against.
data := "data"

# The book the tables are transcribed from. Not in the repository:
# see the extract recipes.
book := "book.pdf"

_default:
    @just --list --unsorted

# ---------------------------------------------------------------- build

# Compile every package.
build:
    go build ./...

# Remove build and coverage artefacts.
clean:
    go clean ./...
    rm -f coverage.out coverage.html

# ----------------------------------------------------------------- test

# Run the whole suite.
test:
    go test ./...

# Run the whole suite, ignoring cached results.
retest:
    go test ./... -count=1

# Run one package's tests verbosely, e.g. `just test-pkg cubic`.
test-pkg pkg:
    go test ./{{ pkg }}/ -v -count=1

# Run tests matching a name, e.g. `just test-run Azeotrope`.
test-run pattern:
    go test ./... -run {{ pattern }} -v -count=1

# Run the suite with the race detector.
race:
    go test ./... -race -count=1

# --------------------------------------------------------------- cover

# Statement coverage for every package.
cover:
    @go test ./... -cover 2>&1 | grep -E '^(ok|---)' || true

# Per-function coverage for one package, e.g. `just cover-pkg cubic`.
cover-pkg pkg:
    go test ./{{ pkg }}/ -coverprofile=coverage.out
    go tool cover -func=coverage.out

# Per-function coverage for one package, only what is under 100%.
cover-gaps pkg:
    @go test ./{{ pkg }}/ -coverprofile=coverage.out > /dev/null
    @go tool cover -func=coverage.out | awk '$3+0 < 100'

# Open a coverage report for one package in the browser.
cover-html pkg:
    go test ./{{ pkg }}/ -coverprofile=coverage.out
    go tool cover -html=coverage.out

# ---------------------------------------------------------------- lint

# Format every package in place.
fmt:
    gofmt -w .

# Report files that are not gofmt-clean, without changing them.
fmt-check:
    @test -z "$(gofmt -l .)" || { echo "not gofmt-clean:"; gofmt -l .; exit 1; }

# Run go vet.
vet:
    go vet ./...

# Everything a change should pass before being committed.
check: fmt-check vet retest examples-build

# ------------------------------------------------------------ generate

# Regenerate every table from the JSON under data/.
generate:
    # The generators are deterministic, so this is a no-op unless the
    # source data has changed.
    go generate ./...

# Show what regenerating would change, without keeping it.
generate-check: generate
    @git diff --stat -- '*/table.go' '*/tables.go' || true

# ------------------------------------------------------------- extract

# Extract book pages into data/, e.g. `just extract 688 689 appendix_c`.
extract start end name:
    # The extracted pages are kept in the repository so that a parsing
    # bug can be traced to its source without locating the book again.
    uv run python -c "import sys; sys.path.insert(0,'scripts'); \
        from extractor import extract_pages; \
        extract_pages('{{ book }}', {{ start }}, {{ end }}, '{{ data }}/{{ name }}.pdf')"

# Re-parse the heat-capacity tables (Appendix C) into data/cp.json.
parse-cp:
    uv run python -c "import sys, json; sys.path.insert(0,'scripts'); \
        from parse import parse_cp_tables; \
        t = parse_cp_tables('{{ data }}/appendix_c_heat_capacities.pdf'); \
        f = open('{{ data }}/cp.json','w',encoding='utf-8'); \
        json.dump(t, f, indent=2); f.write('\n')"

# Re-parse the Antoine constants (Table B.2) into data/b2_antoine.json.
parse-antoine:
    uv run python -c "import sys, json; sys.path.insert(0,'scripts'); \
        from parse import parse_antoine_table; \
        t = parse_antoine_table('{{ data }}/appendix_b2_antoine.pdf'); \
        f = open('{{ data }}/b2_antoine.json','w',encoding='utf-8'); \
        json.dump(t, f, indent=2); f.write('\n')"

# ------------------------------------------------------------ examples

# Build every example.
examples-build:
    go build ./examples/...

# Build and run every example, so a stale one cannot pass unnoticed.
examples:
    #!/usr/bin/env sh
    # The README links to these, and they exercise the public API the way
    # a reader will, which the test suite does not.
    #
    # just only treats a recipe as a script when the shebang is its first
    # line; a comment above it demotes the body to one command per line,
    # where the indentation below is a syntax error.
    set -e
    for dir in examples/*/; do
        name=$(basename "$dir")
        echo ""
        echo "===== $name ====="
        go run "./examples/$name"
    done

# Run one example, e.g. `just example flash`.
example name:
    go run ./examples/{{ name }}

# ----------------------------------------------------------------- doc

# Serve the package documentation locally.
doc:
    @echo "serving on http://localhost:6060/pkg/github.com/rickykimani/zfactor/"
    go run golang.org/x/tools/cmd/godoc@latest -http=:6060

# Print one package's documentation, e.g. `just doc-pkg cubic`.
doc-pkg pkg:
    go doc -all ./{{ pkg }}

# ----------------------------------------------------------------- git

# Summarise how far each branch is from its remote.
status:
    @git status --short
    @git for-each-ref --format='%(refname:short) %(upstream:track)' refs/heads
