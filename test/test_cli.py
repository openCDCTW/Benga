from click.testing import CliRunner
import unittest
import Benga


class CliTest(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

    def test_help(self):
        result = self.runner.invoke(Benga.main, ['-h'])
        self.assertEqual(result.exit_code, 0)

    def tearDown(self):
        self.runner = None


if __name__ == '__main__':
    unittest.main()
