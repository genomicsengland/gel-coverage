import argparse
import json
import logging
import unittest
from gelcoverage.test.output_verifier import OutputVerifier


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser(
        description = 'The JSON verifier checks the output of the GEL'
                      'coverage analyser. It runs format checks and also content'
                      'checks to make sure that the results are coherent.'
                      'This software is intended to be used during the testing'
                      'of the coverage analyser.')
    parser.add_argument('--json', metavar='json',
                        help = 'The JSON file to verify [required]', required = True)
    parser.add_argument('--expected-gene-list', metavar='expected_gene_list',
                        help='The expected comma-separated gene list', default=None)
    args = parser.parse_args()
    coverage_data = json.load(open(args.json, 'r'))

    class TestWrapper(OutputVerifier):
        def run_test(self):
            # It needs to set this config variable so it knows what the settings were
            self.config = coverage_data["parameters"]
            self.verify_output(coverage_data, args.expected_gene_list)

    suite = unittest.TestSuite()
    suite.addTest(TestWrapper("run_test"))
    runner = unittest.TextTestRunner()
    runner.run(suite)


if __name__ == '__main__':
    main()