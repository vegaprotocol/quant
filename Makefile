# Makefile

.PHONY: default
default: test

.PRECIOUS: .testCoverage.html
.testCoverage.html: .testCoverage.txt
	@go tool cover -html="$^" -o "$@"

.PRECIOUS: .testCoverage.txt
.testCoverage.txt:
	@go test -covermode=count -coverprofile="$@" ./...
	@go tool cover -func="$@"

.PHONY: bench
bench: ## Build benchmarking binary (in "$GOPATH/bin"); Run benchmarking
	@go test -run=XXX -bench=. -benchmem -benchtime=1s ./...

.PHONY: build
build:
	@echo "This is a library." ; exit 1

.PHONY: coverage
coverage: .testCoverage.txt ## Generate global code coverage report

.PHONY: coveragehtml
coveragehtml: .testCoverage.html ## Generate global code coverage report in HTML

.PHONY: deps
deps: ## Get the dependencies
	@go mod download

.PHONY: help
help: ## Display this help screen
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

.PHONY: install
install:
	@echo "This is a library." ; exit 1

.PHONY: lint
lint: ## Lint the files
	@go list ./... | xargs -r golint -set_exit_status | sed -e "s#^$$GOPATH/src/##"

.PHONY: msan
msan: ## Run memory sanitizer
	@if ! which clang 1>/dev/null ; then echo "Need clang" ; exit 1 ; fi
	@env CC=clang CGO_ENABLED=1 go test -msan ./...

.PHONY: race
race: ## Run data race detector
	@env CGO_ENABLED=1 go test -race ./...

.PHONY: retest
retest: deps ## Force re-run of all unit tests
	@go test -count=1 ./...

.PHONY: test
test: deps ## Run unit tests
	@go test ./...

.PHONY: vet
vet: ## Run go vet
	@go vet -all ./...

.PHONY: vetshadow
vetshadow: ## Run go vet -shadow
	@go vet -shadow ./...
