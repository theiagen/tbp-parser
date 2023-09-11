from tbp_parser.Parser import Parser
import tbp_parser.tbp_parser as tbp_parser
import pytest
from unittest.mock import patch
import sys

class TestParser:
  @pytest.mark.parametrize(
          'test_args, expected',
          [(['tb_parser.py', '-h'], 0),
          (['tb_parser.py', '--help'], 0),
          (['tb_parser.py', '-v'], 0),
          (['tb_parser.py', '--version'], 0),
          ])

  def test_modules_exit_codes(self, test_args, expected):
      with pytest.raises(SystemExit) as e:
          with patch.object(sys, 'argv', test_args):
              tbp_parser.main()

      assert e.type == SystemExit
      assert e.value.code == expected 