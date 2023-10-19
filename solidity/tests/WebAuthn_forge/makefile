.SILENT:

THIS_FILE := $(lastword $(MAKEFILE_LIST))

define shell-functions
: BEGIN
runcmd() {
	_cmd=$@;

	script_cmd="script -q /dev/null ${_cmd[@]} >&1";
	script -q /dev/null -c echo 2> /dev/null > /dev/null && script_cmd="script -q /dev/null -c \"${_cmd[@]}\" >&1";

	printf "\e[90;1m[\e[90;1mmake: \e[0;90;1mcmd\e[0;90;1m]\e[0m \e[0;93;1m➔ \e[97;1m$_cmd\e[0m\n" \
		&& ( \
			cmd_output=$(eval "$script_cmd" | tee /dev/tty; exit ${PIPESTATUS[0]}); cmd_exit_code=$?; \
			[ -z "$cmd_output" ] || ([ -z "$(tr -d '[:space:]' <<< $cmd_output)" ] && printf "\e[1A"); \
			[[ "$cmd_exit_code" -eq 0 ]] || return $cmd_exit_code \
		) \
		&& printf "\e[032;1m[✔︎] success\e[0m\n\n" \
			|| (_test_exit=$? \
				&& printf "\e[031;1m[✖︎] fail (exit code: $_test_exit)\e[0m\n\n" \
				&& return $_test_exit) \
			&& [ $? -eq 0 ] \
				|| return $?
}
: END
endef

$(shell sed -n '/^: BEGIN/,/^: END/p' $(THIS_FILE) > .make.functions.sh)
SHELL := /bin/bash --init-file .make.functions.sh -i

default:
	printf """\e[37musage:\e[0m\n \
		  \e[90m$$ \e[0;97;1mmake \e[0;92;1mbuild         \e[0;90m➔ \e[32;3mbuild contracts \e[0m\n \
		  \e[90m$$ \e[0;97;1mmake \e[0;92;1mtest          \e[0;90m➔ \e[32;3mrun unit tests \e[0m\n \
		  \e[90m$$ \e[0;97;1mmake \e[0;92;1mlint          \e[0;90m➔ \e[32;3mlint code \e[0m\n \
		  \e[90m$$ \e[0;97;1mmake \e[0;92;1minstall       \e[0;90m➔ \e[32;3minstall dependencies \e[0m\n \
		  \e[90m$$ \e[0;97;1mmake \e[0;92;1mformat        \e[0;90m➔ \e[32;3mapply code style formatting \e[0m\n \
		  \e[90m$$ \e[0;97;1mmake \e[0;92;1mclean         \e[0;90m➔ \e[32;3mremove build artifacts \e[0m\n \
	""" | sed -e 's/^[ \t	]\{1,\}\(.\)/  \1/'

.PHONY: forge-build
forge-build:
	@runcmd forge build

# TODO: differenciate tests w/o ffi and w/ ffi
.PHONY: forge-test
forge-test:
	@runcmd forge test --ffi

# TODO: differenciate tests w/o ffi and w/ ffi
.PHONY: forge-test-verbose
forge-test-verbose:
	@runcmd forge test -vvv --ffi

.PHONY: forge-coverage
forge-coverage:
	@runcmd forge coverage

.PHONY: forge-lint-check
forge-lint-check:
	@runcmd forge fmt test ../../solidity --check

.PHONY: forge-lint
forge-lint:
	@runcmd forge fmt script test ../../solidity

.PHONY: forge-clean
forge-clean:
	@runcmd forge clean

.PHONY: git-submodule-update
git-submodule-update:
	@runcmd git submodule update --init

build: forge-build
contract: build
contracts: build

test: forge-test
tests: test

test-verbose: forge-test-verbose
testv: test-verbose
testsv: test-verbose

coverage: forge-coverage

# output the linting issues
lint: forge-lint-check
linting: lint
linter: lint

# fix the linting issues
lint-write: forge-lint

clean: forge-clean

install: git-submodule-update

#  ℹ️ RUN THIS SCRIPT USING A LEDGER:
#    forge script DeployElliptic.s.sol:MyScript --rpc-url <RPC_URL> --ledger --sender <ACCOUNT_ADDRESS> \
#    [--broadcast]

#    ℹ️ RUN THIS SCRIPT WITH AN ARBITRARY PRIVATE KEY (NOT RECOMMENDED):
#    PRIVATE_KEY=<PRIVATE_KEY> forge script DeployElliptic.s.sol:MyScript --rpc-url <RPC_URL> [--broadcast]
https://ethereum-holesky.publicnode.com

#deploy: forge script DeployElliptic.s.sol:MyScript --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public --private-key b3619b5eb554ae89e6c53ca7d4d1c4dee671a96d16c0a20e597e194609cc2cc2
deploy:
	forge script script/DeployElliptic.s.sol:MyScript --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public --private-key b3619b5eb554ae89e6c53ca7d4d1c4dee671a96d16c0a20e597e194609cc2cc2 --etherscan-api-key HURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5X --verify 

#verify: forge verify-contract --rpc-url https://ethereum-holesky.publicnode.com --etherscan-api-key HURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5XHURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5X 0x1695C55D8bFd9C9035a9197191448381e38d1bF9 DeployElliptic.s.sol:LibraryWrapper 


#forge verify-contract --etherscan-api-key HURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5XHURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5X --chain holesky --force  0x1695C55D8bFd9C9035a9197191448381e38d1bF9 DeployElliptic.s.sol:LibraryWrapper 
forge verify-contract --etherscan-api-key HURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5XHURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5X --chain sepolia --flatten --force  0x3c548295a71383976318819599937fC3daC75289 script/DeployElliptic.s.sol:LibraryWrapper --watch


