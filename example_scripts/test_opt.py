import pytest
pytest.skip(allow_module_level=True)

import sys, os
from pmx import *

options = [
    Option( "-s", "string", "default_string", "Description of what this string options is doing"),
    Option("-f", "float", 3.1415,"Description of what this float option is doing"),
    Option("-i","int",42,"Description of what this int option is doing"),
    Option("-b","bool",True,"Description of what this bool option is doing"),
    Option("-vec1","rvec", [2.5, 3., -1.], "Description of what this real vector option is doing"),
    Option("-vec2","ivec", [2, 3, 0], "Description of what this int vector option is doing"),
    Option("-vec3","svec", ["a", "b", "c"], "Description of what this string vector option is doing")
    ]


file_options = [
    FileOption("-pdb", "r",["pdb"], "protein.pdb", "input pdb file"),
    FileOption("-opdb", "w",["pdb","gro"], "out.pdb", "output pdb or gro file"),
    FileOption("-mpdb", "r/m",["pdb","gro"], "one_of_many.pdb", "several pdb files"),
    ]

help_text = [ "This program does useful things",
              "as long as you use option a,b and c",
              "but not in combination with d and e"]


cmdl = Commandline( sys.argv, options = options, fileoptions = file_options, program_desc = help_text, version = "2.3")



string_opt = cmdl['-s']
input_pdb = cmdl['-pdb']
