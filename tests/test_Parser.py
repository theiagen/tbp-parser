from tbprofiler_parser.Parser import Parser
import tbprofiler_parser.tbprofiler_parser as tbprofiler_parser
import pytest
from unittest.mock import patch
import sys

class TestParser:
  @pytest.mark.parametrize(
          'test_args, expected',
          [(['tbprofiler_parser.py', '-h'], 0),
          (['tbprofiler_parser.py', '--help'], 0),
          (['tbprofiler_parser.py', '-v'], 0),
          (['tbprofiler_parser.py', '--version'], 0),
          ])

  def test_modules_exit_codes(self, test_args, expected):
      with pytest.raises(SystemExit) as e:
          with patch.object(sys, 'argv', test_args):
              tbprofiler_parser.main()

      assert e.type == SystemExit
      assert e.value.code == expected 