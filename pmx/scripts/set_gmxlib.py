#!/usr/bin/env python
import os
from pmx import model


def main():
    path = os.path.abspath(model.__file__)
    dir_path = os.path.dirname(path)
    gmxlib_proteins = os.path.join(dir_path, 'data/mutff45')
    gmxlib_nuc_acids = os.path.join(dir_path, 'data/mutff45dna')

    print('\n  In order to be able to use the hybrid/alchemical force fields \n'
          '  available in pmx, the environment variable GMXLIB needs to be set.\n')

    print('  The path to your pmx force field library for proteins is:')
    print(('  %s\n' % gmxlib_proteins))
    print('  The path to your pmx force field library for nucleic acids is:')
    print(('  %s\n' % gmxlib_nuc_acids))

    print('  Set the relevant GMXLIB path in your shell session as follows:')
    print(('  $ export GMXLIB=%s\n' % gmxlib_proteins))
    print('  Or you can add this directly in your bashrc file.\n')


def entry_point():
    main()


if __name__ == '__main__':
    entry_point()
