import sys
from .view import start

def main(args):
  try:
    str(args[0]) != 'start' and len(args)>1
  except:
    raise ValueError('Invalid command. To start the tool use: python3 -m protein_graph_analyzer start')
  else: start()

if __name__ == '__main__':
    main(sys.argv[1:])
