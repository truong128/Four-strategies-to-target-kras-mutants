"""Generate tamper-evident chemistry file manifests without manufacturing validation.

All validation flags remain false until independently reviewed by qualified chemists.
"""
import argparse
from pathlib import Path
from full.common import digest, jsonout, path_required, GateError


def files_hash_map(filepaths, description):
    result={}
    for file in filepaths:
        source=path_required(file,description)
        if source.name in result:
            raise GateError('Non-unique basename in '+description+': '+source.name)
        result[source.name]=digest(source)
    return result


def prepared(args):
    pdb=path_required(args.pdb,'--pdb')
    output={
        'validated':False,
        'stereochemistry_resolved':False,
        'charge_and_torsion_checked':False,
        'atom_mapping_checked':False,
        'pdb_sha256':digest(pdb),
        'lib_sha256':files_hash_map(args.libs,'--lib'),
        'frcmod_sha256':files_hash_map(args.frcmods,'--frcmod'),
        'independent_verification':'REQUIRED: confirm chemical connectivity, stereochemistry, protonation, parameter benchmarks and alignment of these exact files',
        'warning':'Checksums only prove which bytes are under review. They do NOT prove parameter correctness.'
    }
    jsonout(args.out,output)


def topology(args):
    top=path_required(args.topology,'--topology')
    xyz=path_required(args.coordinates,'--coordinates')
    output={
        'validated':False,
        'stereochemistry_resolved':False,
        'charge_and_torsion_checked':False,
        'atom_mapping_checked':False,
        'topology_sha256':digest(top),
        'coordinates_sha256':digest(xyz),
        'independent_verification':'REQUIRED: review topology atom ordering, water/ions, periodic cell, stereo, charge, equilibration plan and energy check',
        'warning':'This is a checksum scaffold, NOT permission to run modified oligonucleotide MD.'
    }
    jsonout(args.out,output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); sub=parser.add_subparsers(dest='cmd',required=True)
    a=sub.add_parser('parameters');a.add_argument('--pdb',required=True);a.add_argument('--libs',nargs='+',required=True);a.add_argument('--frcmods',nargs='+',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('topology');a.add_argument('--topology',required=True);a.add_argument('--coordinates',required=True);a.add_argument('--out',required=True)
    args=parser.parse_args(); {'parameters':prepared,'topology':topology}[args.cmd](args)
