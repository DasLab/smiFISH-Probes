import sys

print(sys.argv[0])

import sys

try:
    from shlex import quote as cmd_quote
except ImportError:
    from pipes import quote as cmd_quote

cmdline = " ".join(map(cmd_quote, sys.argv[1:]))
print(type(cmdline))