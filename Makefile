.PHONY: test test-unit test-stub test-integration

test: test-unit test-stub

test-unit:
	./tests/test_model_extension_consistency.sh

test-stub:
	./tests/run_stub_tests.sh

test-integration:
	./tests/run_integration_test.sh
