"""Read the supplied Q6 inputs without treating Rosetta/Amber files as EVB terms."""
from pathlib import Path
import gzip
import hashlib
import json
import re

ROOT=Path(__file__).resolve().parents[1]
CASES=('Des27_7_IN','Des27_7_OUT','Des27_WT_IN','Des27_WT_OUT')


def read(path):
    b=Path(path).read_bytes()
    return gzip.decompress(b).decode() if str(path).endswith('.gz') else b.decode()


def sections(text):
    result={};key=None
    for line in text.splitlines():
        line=re.split(r'[#!]',line,maxsplit=1)[0].strip()
        if not line:continue
        if line.startswith('[') and line.endswith(']'):
            key=line[1:-1].lower();result.setdefault(key,[])
        elif key is not None:result[key].append(line.split())
    return result


def pdb_atoms(text):
    atoms={}
    for line in text.splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            serial=int(line[6:11])
            if serial in atoms:raise ValueError('Duplicate PDB atom serial.')
            atoms[serial]={'name':line[12:16].strip(),'residue':line[17:20].strip(),'chain':line[21:22],
                           'resid':line[22:26].strip(),'xyz_angstrom':[float(line[a:a+8]) for a in (30,38,46)]}
    return atoms


def case_paths(case):
    if case not in CASES:raise ValueError(f'Choose one of {CASES}.')
    base=ROOT/'inputs/q6'/case
    return base/(case+'.fep'),base/(case+'.pdb.gz'),base/(case+'.top.gz')


def audit_case(case):
    fp,pp,tp=case_paths(case);fep=sections(read(fp));atoms=pdb_atoms(read(pp));top=read(tp)
    declared=next(line for line in top.splitlines() if 'Total no. of atoms' in line)
    n,nsolute,nwateratoms,_=map(int,declared.split('=')[0].split())
    if len(atoms)!=n:raise ValueError(f'{case}: PDB/topology atom counts disagree.')
    if dict(fep.get('fep',[])).get('states')!='2':raise ValueError('Expected a two-state Q6 model.')
    qm={int(row[0]):int(row[1]) for row in fep['atoms']}
    if len(qm)!=22 or len(set(qm.values()))!=22 or not set(qm.values()).issubset(atoms):raise ValueError('Invalid 22-atom reactive map.')
    charges={int(row[0]):[float(v) for v in row[1:3]] for row in fep['change_charges']}
    if set(charges)!=set(qm):raise ValueError('Every Q atom must have both state charges.')
    totals=[sum(q[i] for q in charges.values()) for i in range(2)]
    if abs(totals[0]-totals[1])>1e-4:raise ValueError('Reactive-region state charges are not conserved.')
    types={int(row[0]):[float(x) for x in row[1:]] for row in fep['bond_types']}
    changes=[]
    for row in fep['change_bonds']:
        i,j,r,p=map(int,row)
        if i not in atoms or j not in atoms:raise ValueError('Bond change references an absent PDB atom.')
        changes.append({'serials':[i,j],'names':[atoms[i]['name'],atoms[j]['name']],
                        'resids':[atoms[i]['resid'],atoms[j]['resid']], 'type_R':r,'type_P':p})
    morse=[]
    for code,values in types.items():
        if len(values)==3:
            d,a,r=values
            if min(values)<=0:raise ValueError('Nonpositive Morse coefficient.')
            morse.append({'type':code,'De_kcal_mol':d,'alpha_per_angstrom':a,'r0_angstrom':r})
    calibrations=[p for p in (ROOT/'inputs').rglob('*') if p.is_file() and 'qfep' in p.name.lower()]
    return {'case':case,'atoms':n,'solute_atoms':nsolute,'water_molecules':(n-nsolute)//nwateratoms,
            'reactive_atoms':len(qm),'reactive_region_charge_R_P_e':totals,'reactive_atom_identities':
            [{'q_index':q,'pdb_serial':serial,**atoms[serial]} for q,serial in qm.items()],
            'morse_terms':morse,'change_bonds':changes,'native_forcefield':'Q-OPLSAA',
            'native_q6_revision':re.search(r'git id=([^\s]+)',top).group(1),
            'calibrated_offset_and_coupling_verified':False,'qfep_analysis_files_present':[str(p.relative_to(ROOT)) for p in calibrations],
            'scope':'Static published-input audit, not energy parity or runnable full EVB approval.',
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (fp,pp,tp)}}


def audit_all():
    manifest=json.loads((ROOT/'inputs/source_manifest.json').read_text())
    for item in manifest['files']:
        if hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()!=item['stored_sha256']:
            raise ValueError(f'Input hash differs: {item["path"]}')
    return {'schema':'kemp_published_input_audit_v1','cases':[audit_case(c) for c in CASES],
            'empty_source_files':manifest['empty_archive_members_omitted'],
            'full_paper_evb_reproduction_ready':False,
            'blockers':['Empty source window/equilibration inputs and generic path/restart placeholders.',
                        'Calibrated state offset/H12 and Qfep reconstruction settings not verified in the supplied archive.',
                        'Q6/OPLS SCAAS/local-reaction-field to OpenMM whole-Hamiltonian parity not established.'],
            'live_openmm_modes':['Amber-input compatibility smoke (classical, vacuum restrained).','Published Q6 Morse component energy/force replay (not a whole-enzyme PMF).']}
