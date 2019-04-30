# Makefile

.PHONY: default bench built coverage coveragehtml deps help install lint msan race test vet

default: test

.PRECIOUS: .testCoverage.html
.testCoverage.html: .testCoverage.txt
	@go tool cover -html="$^" -o "$@"

.PRECIOUS: .testCoverage.txt
.testCoverage.txt:
	@go test -covermode=count -coverprofile="$@" ./...
	@go tool cover -func="$@"

bench: ## Build benchmarking binary (in "$GOPATH/bin"); Run benchmarking
	@go test -run=XXX -bench=. -benchmem -benchtime=1s ./...

build:
	@echo "This is a library." ; exit 1

coverage: .testCoverage.txt ## Generate global code coverage report

coveragehtml: .testCoverage.html ## Generate global code coverage report in HTML

deps: ## Get the dependencies
	@go mod download

help: ## Display this help screen
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

install:
	@echo "This is a library." ; exit 1

lint: ## Lint the files
	@go list ./... | xargs -r golint -set_exit_status | sed -e "s#^$$GOPATH/src/##"

msan: ## Run memory sanitizer
	@if ! which clang 1>/dev/null ; then echo "Need clang" ; exit 1 ; fi
	@env CC=clang CGO_ENABLED=1 go test -msan ./...

race: ## Run data race detector
	@env CGO_ENABLED=1 go test -race ./...

test: deps ## Run unit tests
	@go test ./...

vet: ## Run go vet
	@go vet -all ./...
