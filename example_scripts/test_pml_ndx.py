import sys, os
from pmx import *
from pmx.ndx import *

import pytest
pytest.skip(allow_module_level=True)
# skip because pymol not installed as dependency

from pymol import cmd, stored


def write_ndx( sel="(all)", fname = "index.ndx" ):
    name = sel
    stored.ID = []
    cmd.iterate(sel,"stored.ID.append(ID)")
    ndx_group = IndexGroup( name, ids = stored.ID )
    if os.path.isfile( fname ):
        # append to existing file
        ndx_file = IndexFile(fname)
        ndx_file.add_group( ndx_group )
    else:
        ndx_file = IndexFile( groups = [ndx_group] ) # new
    ndx_file.write( fname )

def __sel_from_id_list( name, ids ):
    if ids:
        cmd.select( name, "ID %d" % ids[0] )
        for idx in ids[1:]:
            cmd.select( name, "%s or ID %d" % (name, idx) )

def load_ndx( fname = "index.ndx", names = []):
    if not os.path.isfile( fname ): return
    ndx_file = IndexFile( fname )
    if not names:
        names = ndx_file.names
    for name in names:
        print(name)
        if name in ndx_file.dic:
            ids = ndx_file[name].ids
            __sel_from_id_list( name, ids )

    cmd.deselect()

cmd.extend("write_ndx",write_ndx)
cmd.extend("load_ndx",load_ndx)
