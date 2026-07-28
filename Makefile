.PHONY: test test-stub test-integration

test: test-stub

test-stub:
	./tests/run_stub_tests.sh

test-integration:
	./tests/run_integration_test.sh
