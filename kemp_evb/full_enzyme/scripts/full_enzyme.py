"""Full native-input EVBOpenMM adaptation; finite restrained droplet, not Q6 LRF.

No fitted coupling/offset or physical barrier is assumed. The validation command
checks construction and finite energies; it does not establish force-field parity.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import app, unit
from kemp_inputs import ROOT, CASES, case_paths, read, sections
from q_topology import parse_topology


def canonical(atoms):
    t=tuple(atoms)
    return min(t,t[::-1])


def patched_terms(data, fep, kind, state):
    width={'bonds':2,'angles':3,'torsions':4,'impropers':4}[kind]
    singular={'bonds':'bond','angles':'angle','torsions':'torsion','impropers':'improper'}[kind]
    changes=fep['change_'+kind]
    replaced={canonical(map(int,row[:width])) for row in changes}
    result=[]
    for row in data['bonded'][kind]['terms']:
        if canonical(row[:-1]) not in replaced:
            result.append((row[:-1],data['bonded'][kind]['parameter_types'][row[-1]],False))
    types={int(row[0]):list(map(float,row[1:])) for row in fep[singular+'_types']}
    for row in changes:
        code=int(row[width+state])
        if code:result.append((list(map(int,row[:width])),types[code],True))
    return result


def state_system(data,fep,state):
    """Parameter carrier: geometric LJ is constructed by EVBOpenMM, not this NB force."""
    system=mm.System();n=data['atoms']
    for i in data['atom_type_codes']:system.addParticle(data['type_masses_da'][i-1])
    bonds=mm.HarmonicBondForce();angles=mm.HarmonicAngleForce();torsions=mm.PeriodicTorsionForce()
    morse={};harmonic_impropers=[];edges=set()
    for atoms,p,changed in patched_terms(data,fep,'bonds',state):
        i,j=[a-1 for a in atoms];edges.add(tuple(sorted((i,j))))
        if changed and len(p)==3:
            de,alpha,r0=p;morse[(i,j)]=(de,alpha,r0)
            # Curvature carrier, replaced with Morse in the EVB scaffold.
            k=2*de*alpha**2
        else:k,r0=p[:2]
        bonds.addBond(i,j,r0*.1,k*418.4)
    for atoms,p,changed in patched_terms(data,fep,'angles',state):
        angles.addAngle(*[a-1 for a in atoms],math.radians(p[1]),p[0]*4.184)
    for kind in ('torsions','impropers'):
        for atoms,p,changed in patched_terms(data,fep,kind,state):
            ix=[a-1 for a in atoms]
            if kind=='impropers' and changed:
                harmonic_impropers.append((ix,p));continue
            if kind=='impropers':k,phase=p[:2];period=2
            else:
                k,period,phase=p[:3]
                if not changed:k/=p[3]
            torsions.addTorsion(*ix,int(abs(period)),math.radians(phase),k*4.184)
    for force in (bonds,angles,torsions):system.addForce(force)
    nb=mm.NonbondedForce();nb.setNonbondedMethod(nb.NoCutoff)
    q=list(data['charges_e']);coeff=[]
    tables=data['native_LJ_tables']
    for code in data['atom_type_codes']:
        coeff.append([tables[key][code-1] for key in ('sqrt (Aii) normal:','sqrt (Bii) normal:','sqrt (Aii) 1-4:','sqrt (Bii) 1-4:')])
    qmap={int(row[0]):int(row[1])-1 for row in fep['atoms']}
    types={row[0]:list(map(float,row[1:])) for row in fep['atom_types']}
    for row in fep['change_atoms']:
        p=types[row[state+1]];coeff[qmap[int(row[0])]]=[p[0],p[1],p[4],p[5]]
    for row in fep['change_charges']:q[qmap[int(row[0])]]=float(row[state+1])
    def lj(a,b):
        if a==0 or b==0:return 1.,0.
        return (a/b)**(1/3)*.1,(b**4/(4*a*a))*4.184
    for charge,c in zip(q,coeff):nb.addParticle(charge,*lj(c[0],c[1]))
    # Preserve compact native host exceptions; rebuild all exceptions involving
    # Q atoms from each state's graph to follow bond/angle/1-4 changes.
    excluded={tuple(sorted((i-1,j-1))) for key in ('short_pairs','long_pairs') for i,j in data['pairs']['excluded'][key]}
    pairs14={tuple(sorted((i-1,j-1))) for key in ('short_pairs','long_pairs') for i,j in data['pairs']['one_four'][key]}
    active=set(qmap.values())
    excluded={p for p in excluded if not active.intersection(p)}
    pairs14={p for p in pairs14 if not active.intersection(p)}
    graph=[set() for _ in range(n)]
    for i,j in edges:graph[i].add(j);graph[j].add(i)
    for i in active:
        seen={i};front={i}
        for depth in range(1,4):
            front=set().union(*(graph[j] for j in front))-seen if front else set();seen|=front
            for j in front:
                pair=tuple(sorted((i,j)))
                (excluded if depth<3 else pairs14).add(pair)
    pairs14-=excluded
    for i,j in sorted(excluded):nb.addException(i,j,0.,1.,0.)
    for i,j in sorted(pairs14):
        a=coeff[i][2]*coeff[j][2];b=coeff[i][3]*coeff[j][3]
        sigma=(a/b)**(1/6)*.1 if a and b else 1.
        epsilon=b*b/(4*a)*4.184 if a and b else 0.
        nb.addException(i,j,q[i]*q[j]*data['coulomb_14_scale'],sigma,epsilon)
    system.addForce(nb)
    # Native rigid three-site water; no reactive constraints.
    for o in range(data['solute_atoms'],n,3):
        system.addConstraint(o,o+1,.0957);system.addConstraint(o,o+2,.0957);system.addConstraint(o+1,o+2,.15183)
    return system,morse,harmonic_impropers,edges


def build(case,out):
    from evbopenmm import EVBSystem
    from evbopenmm.models import MorseBondParameter
    fp,pp,tp=case_paths(case);data=parse_topology(read(tp));fep=sections(read(fp))
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    pdb_path=out/'native_identity.pdb'
    # Q uses ASP for neutral and charged Asp; OpenMM requires distinct names.
    native=read(pp);rows=native.splitlines();neutral_asp={line[21:26] for line in rows if line.startswith(('ATOM  ','HETATM')) and line[17:20]=='ASP' and line[12:16].strip() in ('HD1','HD2')}
    aliases={'ASH':'DPN','GLH':'EPN','HID':'H1N','HIE':'H2N','HIP':'H3N'}
    rows=[line[:17]+aliases.get(line[17:20],line[17:20])+line[20:] if line.startswith(('ATOM  ','HETATM')) else line for line in rows]
    rows=[line[:17]+'DPN'+line[20:] if line.startswith(('ATOM  ','HETATM')) and line[17:20]=='ASP' and line[21:26] in neutral_asp else line for line in rows]
    rows=[line[:21]+'A'+line[22:] if line.startswith(('ATOM  ','HETATM')) else line for line in rows]
    pdb_path.write_text('\n'.join(rows)+'\n')
    systems={};morse={};impropers={};graphs={}
    for state,label in enumerate(('R','P')):
        systems[label],morse[label],impropers[label],graphs[label]=state_system(data,fep,state)
    pdb=app.PDBFile(str(pdb_path));atoms=list(pdb.topology.atoms())
    def ref(label,i):
        a=atoms[i];return (label,a.residue.chain.id,int(a.residue.id),a.name)
    mapping={ref('R',i):ref('P',i) for i in range(len(atoms))}
    if len(mapping)!=data['atoms']:raise ValueError('Native PDB has ambiguous atom identities.')
    evb=EVBSystem({'R':pdb_path,'P':pdb_path},force_field='opls-aam',lj_combining='geometric',atom_mapping=mapping,
                  prebuilt_state_systems=systems,auto_build_forces=False,reactive_residue_names={'SUB','ASP'})
    qindices=[int(row[1])-1 for row in fep['atoms']]
    evb.set_qm_atoms([ref('R',i) for i in qindices])
    params={}
    for state in ('R','P'):
        for pair,(de,alpha,r0) in morse[state].items():
            params[tuple(ref('R',i) for i in pair)]=MorseBondParameter(r0*.1,de*4.184,alpha*10)
    detected=evb.get_forming_breaking_bond_pairs()
    normalized={evb._normalize_reactive_pair(k):v for k,v in params.items()}
    if set(detected)!=set(normalized):raise ValueError('Reactive bond detection differs from native Morse pairs.')
    evb.build_reactive_forces(reactive_morse_parameters=normalized,
        soft_pair_parameters_by_atom_type={'H':(5.,2.5),'C':(90.,2.5),'N':(90.,2.5),'O':(90.,2.5)},
        state_dependent_exclusion_patches=True,mm_cutoff_distance=9.9,qm_cutoff_distance=9.9)
    # PDB standard-residue guessing loses native protonation-variant bonds.
    # Rebuild topology connectivity from each imported state's actual bond graph.
    for label in ('R','P'):
        topology=evb.states[label].topology;local=list(topology.atoms());topology._bonds=[]
        for i,j in sorted(graphs[label]):
            # Native zero-k H-H water entries carry a SHAKE distance, not a chemical bond.
            if i>=data['solute_atoms'] and (i-data['solute_atoms'])%3!=0 and (j-data['solute_atoms'])%3!=0:continue
            topology.addBond(local[i],local[j])
    # Changed Q improper has a harmonic wrapped form, not the host periodic form.
    for label in ('R','P'):
        if not impropers[label]:continue
        parameter=evb._state_lambda_parameter_names[label]
        force=mm.CustomTorsionForce(f'{parameter}*0.5*k*delta^2;delta=atan2(sin(theta-theta0),cos(theta-theta0))')
        defaults=[f.getGlobalParameterDefaultValue(j) for f in evb.evb_scaffold.getForces() if hasattr(f,'getNumGlobalParameters') for j in range(f.getNumGlobalParameters()) if f.getGlobalParameterName(j)==parameter]
        force.addGlobalParameter(parameter,defaults[0])
        force.addPerTorsionParameter('k');force.addPerTorsionParameter('theta0')
        for ix,p in impropers[label]:force.addTorsion(*ix,[p[0]*4.184,math.radians(p[1])])
        evb.evb_scaffold.addForce(force);evb._register_force_to_group(force,term='native_harmonic_improper',state=label)
    # Explicit changed boundary: finite droplet with flat-bottom oxygen wall.
    cx,cy,cz=[v*.1 for v in data['solvent_center_angstrom']]
    wall=mm.CustomExternalForce('0.5*k*max(0,sqrt((x-cx)^2+(y-cy)^2+(z-cz)^2)-radius)^2')
    for key,value in {'k':10000.,'cx':cx,'cy':cy,'cz':cz,'radius':data['solvent_radius_angstrom']*.1}.items():wall.addGlobalParameter(key,value)
    for i in range(data['solute_atoms'],data['atoms'],3):wall.addParticle(i,[])
    evb.evb_scaffold.addForce(wall);evb._register_force_to_group(wall,term='droplet_wall',state=None)
    for i in data['excluded_atom_indices_1based']:evb.evb_scaffold.setParticleMass(i-1,0.)
    # The library's zero-min Morse convention changes P-R by -44.5 kcal/mol.
    # This correction is for energy bookkeeping only, not a fitted state offset.
    correction=(sum(v[0] for v in morse['R'].values())-sum(v[0] for v in morse['P'].values()))
    metadata={'case':case,'particles':data['atoms'],'reactive_atoms':len(qindices),'reactive_pairs':len(detected),
       'boundary':'finite nonperiodic droplet; direct Coulomb within 9.9 nm; harmonic water-O wall; native outer protein frozen',
       'native_q6_reproduction':False,'calibrated':False,'morse_P_minus_R_correction_kcal_mol':correction,
       'limitations':['No LRF/SCAAS equivalence.','No validated coupling or physical state offset.','No converged free-energy barrier.','State exception conversion requires separate numerical validation.']}
    return evb,data,metadata


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=CASES);p.add_argument('--output-dir',type=Path)
    p.add_argument('--mode',choices=('validate','trial','prepare','relax','equilibrate','production','run','diagnose','profile','target-fit','scan'),default='validate')
    p.add_argument('--config',type=Path);p.add_argument('--resume',action='store_true')
    p.add_argument('--records',type=Path,help='Saved full-enzyme raw JSON for profile/diagnose only; never resumes MD.')
    p.add_argument('--max-runtime-seconds',type=float)
    p.add_argument('--coupling',type=float);p.add_argument('--offset',type=float,default=0.)
    p.add_argument('--bins',type=int);p.add_argument('--min-points',type=int,default=20)
    p.add_argument('--target-barrier',type=float)
    p.add_argument('--skip-samples',type=int,default=0)
    p.add_argument('--open',action='store_true')
    p.add_argument('--couplings',type=float,nargs='+',default=[60.,80.,100.4,120.,140.])
    p.add_argument('--offsets',type=float,nargs='+',default=[-130.,-120.,-110.,-100.,-90.])
    p.add_argument('--no-refine',action='store_true')
    a=p.parse_args()
    if a.records is not None and a.mode not in ('profile','diagnose'):p.error('--records is for profile/diagnose only.')
    if a.mode in ('run','prepare','relax','equilibrate','production') and a.coupling is not None:p.error('--coupling is a profile-analysis option; mapping sampling saves both diabatic energies.')
    if a.bins is None:a.bins=18 if a.mode=='scan' else 40
    if a.max_runtime_seconds is not None and (not math.isfinite(a.max_runtime_seconds) or a.max_runtime_seconds<=0):p.error('--max-runtime-seconds must be positive and finite')
    if a.mode not in ('validate','trial'):
        from full_workflow import run
        return run(a)
    a.case=a.case or CASES[0];a.output_dir=a.output_dir or ROOT/'runs/full_enzyme_validation'
    evb,data,meta=build(a.case,a.output_dir)
    integrator=mm.LangevinMiddleIntegrator(300.,1.,.0001);integrator.setRandomNumberSeed(101)
    platform=mm.Platform.getPlatformByName('CPU')
    simulation=app.Simulation(evb.states['R'].topology,evb.evb_scaffold,integrator,platform,{'Threads':'2'})
    simulation.context.setPositions(np.array(data['coordinates_angstrom'])*.1)
    checks=[]
    for lam in (1.,.5,0.):
        evb.set_state_lambdas(simulation,[lam])
        state=simulation.context.getState(getEnergy=True,getForces=True)
        energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces=state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
        if not math.isfinite(energy) or not np.isfinite(forces).all():raise ValueError('Nonfinite full-system energy/forces.')
        checks.append({'lambda_R':lam,'potential_energy_kj_mol':energy,'max_force_kj_mol_nm':float(np.linalg.norm(forces,axis=1).max())})
    # One-coordinate force/energy consistency check (CPU precision tolerance).
    evb.set_state_lambdas(simulation,[.5])
    xyz=np.array(data['coordinates_angstrom'])*.1
    particle=int(sections(read(case_paths(a.case)[0]))['atoms'][0][1])-1
    simulation.context.setPositions(xyz)
    analytic=simulation.context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)[particle,0]
    h=.0001;energies=[]
    for direction in (-1,1):
        shifted=xyz.copy();shifted[particle,0]+=direction*h;simulation.context.setPositions(shifted)
        energies.append(simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    numerical=-(energies[1]-energies[0])/(2*h)
    error=abs(analytic-numerical)
    if error>max(5.,abs(analytic)*.02):raise ValueError(f'Force finite-difference check failed: {analytic} versus {numerical}')
    simulation.context.setPositions(xyz)
    meta['force_finite_difference_check']={'particle_zero_based':particle,'component':'x','analytic_kj_mol_nm':float(analytic),'numerical_kj_mol_nm':float(numerical),'absolute_error_kj_mol_nm':float(error),'scope':'one coordinate; not independent native Hamiltonian parity'}
    meta['input_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in case_paths(a.case)}
    meta['driver_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    meta['checks']=checks
    residual=checks[1]['potential_energy_kj_mol']-.5*(checks[0]['potential_energy_kj_mol']+checks[2]['potential_energy_kj_mol'])
    if abs(residual)>.03:raise ValueError(f'Mapping linearity failed: {residual} kJ/mol')
    meta['mapping_linearity_residual_kj_mol']=residual
    meta['status']='finite construction check; independent Hamiltonian parity and calibration still required'
    if a.mode=='trial':
        evb.set_state_lambdas(simulation,[1.])
        simulation.context.setVelocitiesToTemperature(300.,101)
        with (a.output_dir/'topology.pdb').open('w') as handle:
            app.PDBFile.writeFile(simulation.topology,simulation.context.getState(getPositions=True).getPositions(),handle)
        simulation.reporters.append(app.DCDReporter(str(a.output_dir/'trajectory.dcd'),10))
        records=[]
        for lam in (1.,.5,0.):
            evb.set_state_lambdas(simulation,[lam]);simulation.step(10)
            for sample in range(2):
                simulation.step(10)
                energies=evb._sample_diabatic_state_energies_kj_mol(simulation)
                values={key:float(value) for key,value in zip(('R','P'),energies)}
                if not all(math.isfinite(v) for v in values.values()):raise ValueError('Nonfinite trial diabat energy.')
                records.append({'lambda_R':lam,'step':simulation.currentStep,'diabatic_energies_kj_mol_zero_min_morse':values})
        raw={'schema':'kemp_full_enzyme_engineering_trial_v1','metadata':meta,'samples':records,
             'calibrated':False,'physical_barrier_available':False,'production_steps_per_window':20,'equilibration_steps_per_window':10,'timestep_fs':.1}
        (a.output_dir/'replica_101_raw.json').write_text(json.dumps(raw,indent=2,allow_nan=False)+'\n')
        with (a.output_dir/'energies.csv').open('w') as handle:
            writer=csv.writer(handle);writer.writerow(['step','lambda_R','E_R_kj_mol','E_P_kj_mol'])
            for row in records:writer.writerow([row['step'],row['lambda_R'],row['diabatic_energies_kj_mol_zero_min_morse']['R'],row['diabatic_energies_kj_mol_zero_min_morse']['P']])
        meta['trial']={'total_steps':simulation.currentStep,'samples':len(records),'status':'engineering only; insufficient for free-energy inference'}
    (a.output_dir/'validation.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
    print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
