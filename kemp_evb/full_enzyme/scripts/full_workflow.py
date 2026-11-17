"""Restartable relaxation, thermal equilibration and mapping production.

Raw mapping production does not require a coupling. This is an explicitly
changed, uncalibrated Hamiltonian, not a reproduced Q6 activation profile.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time
import uuid
import numpy as np
import openmm as mm
from openmm import app, unit
from kemp_inputs import ROOT, CASES, case_paths
from full_enzyme import build

DEFAULTS=dict(case=CASES[0],platform='CPU',cpu_threads=2,seed=101,
    temperature_kelvin=300.,max_temperature_kelvin=1500.,timestep_fs=.25,friction_per_ps=1.,
    relaxation_steps=100,relaxation_lbfgs_iterations=200,relaxation_max_displacement_nm=.0005,
    relaxation_force_tolerance_kj_mol_nm=250.,
    heating_temperatures_kelvin=[50.,150.,300.],heating_steps_per_temperature=200,
    equilibration_steps=1000,lambda_schedule=[float(x) for x in np.linspace(1,0,21)],
    window_equilibration_steps=200,production_steps_per_window=1000,
    sample_interval=20,trajectory_interval=100,checkpoint_interval=100)


def load_config(path):
    supplied=json.loads(Path(path).read_text()) if path else {}
    unknown=set(supplied)-set(DEFAULTS)-{'schema','note'}
    if unknown:raise ValueError(f'Unknown configuration keys: {sorted(unknown)}')
    cfg={**DEFAULTS,**supplied}
    if cfg.get('schema','kemp_full_mapping_v1')!='kemp_full_mapping_v1':raise ValueError('Unsupported config schema.')
    if cfg['case'] not in CASES:raise ValueError('Unknown native case.')
    for key in ('cpu_threads','seed','relaxation_steps','relaxation_lbfgs_iterations','heating_steps_per_temperature','equilibration_steps','window_equilibration_steps','production_steps_per_window','sample_interval','trajectory_interval','checkpoint_interval'):
        value=cfg[key]
        if isinstance(value,bool) or not isinstance(value,int) or value<=0:raise ValueError(f'{key} must be a positive integer.')
    for key in ('temperature_kelvin','max_temperature_kelvin','timestep_fs','friction_per_ps','relaxation_max_displacement_nm','relaxation_force_tolerance_kj_mol_nm'):
        value=cfg[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:raise ValueError(f'{key} must be finite and positive.')
    if cfg['max_temperature_kelvin']<=cfg['temperature_kelvin']:raise ValueError('max_temperature_kelvin must exceed temperature_kelvin.')
    if cfg['timestep_fs']>1.:raise ValueError('Flexible reactive bonds require timestep_fs <= 1; validate stability before increasing.')
    if cfg['relaxation_max_displacement_nm']>.002:raise ValueError('Use relaxation_max_displacement_nm <= 0.002 for bounded minimization.')
    schedule=cfg['lambda_schedule']
    if not isinstance(schedule,list) or len(schedule)<2 or any(isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x) or x<0 or x>1 for x in schedule):raise ValueError('lambda_schedule must contain finite weights from 1 to 0.')
    if schedule[0]!=1 or schedule[-1]!=0 or any(a<=b for a,b in zip(schedule,schedule[1:])):raise ValueError('lambda_schedule must decrease strictly from 1 to 0.')
    temps=cfg['heating_temperatures_kelvin']
    if not isinstance(temps,list) or not temps or any(isinstance(t,bool) or not isinstance(t,(float,int)) or not math.isfinite(t) or t<=0 for t in temps) or temps[-1]!=cfg['temperature_kelvin']:raise ValueError('Heating temperatures must be positive and end at temperature_kelvin.')
    if cfg['checkpoint_interval']%cfg['sample_interval']:raise ValueError('checkpoint_interval must be a multiple of sample_interval.')
    if cfg['sample_interval']>cfg['production_steps_per_window']:raise ValueError('sample_interval exceeds production_steps_per_window.')
    return cfg


def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)


def fingerprint(cfg):
    paths=[Path(__file__),Path(__file__).with_name('full_enzyme.py'),Path(__file__).with_name('q_topology.py'),*case_paths(cfg['case'])]
    from evbopenmm import system
    paths.append(Path(system.__file__))
    payload={'config':cfg,'sources':{str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'openmm_version':mm.version.version}
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest(),payload


def finite_state(sim):
    state=sim.context.getState(getEnergy=True,getForces=True,getPositions=True,getVelocities=True)
    energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    kinetic=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    positions=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    forces=state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
    velocities=state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
    if not math.isfinite(energy) or not math.isfinite(kinetic) or not all(np.isfinite(a).all() for a in (positions,forces,velocities)):raise ValueError('Nonfinite energy, forces, positions or velocities; stage stopped.')
    return state,energy,positions,forces


def bounded_minimize(sim,cfg,output,stop=None):
    """Projected steepest descent with Armijo-like energy decrease and capped motion.

    Each line-search candidate is constrained, checked, and either accepted or
    rolled back. No invalid iterate is passed to heating. The iteration cap is
    reported separately from reaching the force tolerance.
    """
    mobile=np.array([sim.system.getParticleMass(i).value_in_unit(unit.dalton)>0 for i in range(sim.system.getNumParticles())])
    sim.context.applyConstraints(1e-6)
    trace=[];reason='iteration_limit';start_energy=None
    for iteration in range(cfg['relaxation_steps']+1):
        state,energy,x,f=finite_state(sim)
        if start_energy is None:start_energy=energy
        f[~mobile]=0;fmax=float(np.linalg.norm(f,axis=1).max())
        trace.append({'iteration':iteration,'potential_energy_kj_mol':energy,'max_mobile_force_kj_mol_nm':fmax})
        if iteration%20==0:print(f'Relax {iteration}/{cfg["relaxation_steps"]}: E={energy:.3f} kJ/mol, max mobile force={fmax:.2f}',flush=True)
        if fmax<=cfg['relaxation_force_tolerance_kj_mol_nm']:reason='force_tolerance';break
        if iteration==cfg['relaxation_steps']:break
        if stop and stop():reason='time_limit';break
        displacement=cfg['relaxation_max_displacement_nm']*f/max(fmax,1.)
        accepted=False
        for _ in range(24):
            try:
                sim.context.setPositions(x+displacement);sim.context.applyConstraints(1e-6)
                candidate,candidate_energy,_,_=finite_state(sim)
                if candidate_energy<energy:
                    accepted=True;break
            except (ValueError,mm.OpenMMException):pass
            sim.context.setPositions(x);displacement*=.5
        if not accepted:
            sim.context.setPositions(x);reason='line_search_stalled';break
    descent_iterations=trace[-1]['iteration'];lbfgs_attempted=False;lbfgs_calls=0
    if reason!='time_limit' and reason!='force_tolerance' and cfg.get('relaxation_lbfgs_iterations',0):
        lbfgs_attempted=True;backup=finite_state(sim)[2].copy();before=finite_state(sim)[1]
        class DeadlineReporter(mm.MinimizationReporter):
            def report(self,iteration,x,grad,arguments):
                self.calls+=1
                if iteration%50==0:print(f'L-BFGS iteration {iteration}',flush=True)
                return bool(stop and stop())
        try:
            reporter=DeadlineReporter();reporter.calls=0
            sim.minimizeEnergy(tolerance=cfg['relaxation_force_tolerance_kj_mol_nm'],maxIterations=cfg['relaxation_lbfgs_iterations'],reporter=reporter)
            lbfgs_calls=reporter.calls
            _,energy,_,forces=finite_state(sim);forces[~mobile]=0
            if energy>before+.05:raise ValueError('L-BFGS increased constrained potential energy.')
            fmax=float(np.linalg.norm(forces,axis=1).max())
            trace.append({'iteration':trace[-1]['iteration']+1,'potential_energy_kj_mol':energy,'max_mobile_force_kj_mol_nm':fmax})
            reason='time_limit' if stop and stop() else ('force_tolerance' if fmax<=cfg['relaxation_force_tolerance_kj_mol_nm'] else 'bounded_lbfgs_completed')
        except (ValueError,mm.OpenMMException):
            sim.context.setPositions(backup)
            raise ValueError('L-BFGS failed after bounded warm-up; heating was not started.')
    result={'algorithm':'displacement-limited constrained descent followed by L-BFGS' if lbfgs_attempted else 'displacement-limited constrained steepest descent with energy-decreasing line search',
        'stop_reason':reason,'converged_force_tolerance':reason=='force_tolerance','initial_energy_kj_mol':start_energy,
        'final_energy_kj_mol':trace[-1]['potential_energy_kj_mol'],'descent_iterations':descent_iterations,'lbfgs_iteration_limit':cfg.get('relaxation_lbfgs_iterations',0),'lbfgs_report_calls':lbfgs_calls,
        'final_max_mobile_force_kj_mol_nm':trace[-1]['max_mobile_force_kj_mol_nm'],
        'force_note':'Includes constraint-normal forces; bounded relaxation is not proof of a stationary unconstrained minimum.'}
    atomic_json(output/'relaxation.json',result)
    with (output/'relaxation.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(trace[0]));writer.writeheader();writer.writerows(trace)
    if reason=='line_search_stalled' and result['descent_iterations']==0:raise ValueError('Relaxation made no progress; inspect relaxation.json.')
    return result


class Store:
    """Samples and their exact restart cursor/checkpoint commit in one transaction."""
    def __init__(self,path,signature,resume):
        existed=Path(path).exists()
        if existed and not resume:raise ValueError('Run already exists. Use --resume for identical settings or choose a new --output-dir.')
        self.db=sqlite3.connect(path)
        self.db.execute('CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value BLOB)')
        self.db.execute('CREATE TABLE IF NOT EXISTS samples (window INTEGER, local_step INTEGER, step INTEGER, lambda_R REAL, E_R REAL, E_P REAL, mapping_energy REAL, PRIMARY KEY(window,local_step))')
        previous=self.get('fingerprint')
        if previous is not None and previous.decode()!=signature:
            self.db.close()
            raise ValueError('Config, inputs, OpenMM or force code changed. Use a new output directory; checkpoints cannot be mixed.')
        self.db.execute('INSERT OR IGNORE INTO meta VALUES (?,?)',('fingerprint',signature.encode()));self.db.commit()
    def get(self,key):
        row=self.db.execute('SELECT value FROM meta WHERE key=?',(key,)).fetchone();return row[0] if row else None
    def cursor(self):
        value=self.get('cursor');return json.loads(value) if value else {'phase':'prepared'}
    def save(self,sim,cursor,samples=()):
        checkpoint=sim.context.createCheckpoint();xml=mm.XmlSerializer.serialize(sim.context.getState(getPositions=True,getVelocities=True,getParameters=True,getIntegratorParameters=True)).encode()
        with self.db:
            self.db.executemany('INSERT INTO samples VALUES (?,?,?,?,?,?,?)',samples)
            for key,value in (('cursor',json.dumps(cursor).encode()),('checkpoint',checkpoint),('state_xml',xml)):
                self.db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)',(key,value))
    def restore(self,sim):
        if self.get('checkpoint'):
            try:sim.context.loadCheckpoint(self.get('checkpoint'))
            except mm.OpenMMException as exc:raise ValueError('Checkpoint requires the same OpenMM platform/hardware; use a new run for another platform.') from exc
            sim.currentStep=int(round(sim.context.getStepCount()))
    def export(self,out,cfg,meta):
        rows=self.db.execute('SELECT window,local_step,step,lambda_R,E_R,E_P,mapping_energy FROM samples ORDER BY window,local_step').fetchall()
        with (out/'energies.csv').open('w',newline='') as handle:
            writer=csv.writer(handle);writer.writerow(['window','window_production_step','step','lambda_R','E_R_kj_mol','E_P_kj_mol','mapping_energy_kj_mol']);writer.writerows(rows)
        windows=[]
        for i,lam in enumerate(cfg['lambda_schedule']):
            selected=[row for row in rows if row[0]==i]
            windows.append({'index':i,'lambda_R':lam,'samples':[{'step':row[2],'E_R_kj_mol':row[4],'E_P_kj_mol':row[5]} for row in selected]})
        atomic_json(out/f'replica_{cfg["seed"]}_raw.json',{'schema':'kemp_full_mapping_records_v1','metadata':meta,'config':cfg,'cursor':self.cursor(),
            'morse_energy_convention':'EVBOpenMM zero-minimum; add +44.5 kcal/mol to P relative to R for the native convention only',
            'calibrated':False,'physical_barrier_available':False,'sample_count':len(rows),'windows':windows})
        atomic_json(out/'progress.json',self.cursor())
        if self.get('state_xml'):(out/'latest_state.xml').write_bytes(self.get('state_xml'))


def diagnose(out,seed,records=None):
    raw=json.loads((Path(records) if records else out/f'replica_{seed}_raw.json').read_text())
    temperature=raw['config']['temperature_kelvin'];beta=1/(.00831446261815324*temperature)
    windows=raw['windows'];edges=[]
    def logmeanexp(values):
        values=np.asarray(values);maximum=values.max();return float(maximum+np.log(np.exp(values-maximum).mean()))
    for left,right in zip(windows,windows[1:]):
        a=np.array([row['E_R_kj_mol']-row['E_P_kj_mol'] for row in left['samples']])
        b=np.array([row['E_R_kj_mol']-row['E_P_kj_mol'] for row in right['samples']])
        edge={'left_window':left['index'],'right_window':right['index'],'sample_counts':[len(a),len(b)]}
        if len(a) and len(b):
            delta=right['lambda_R']-left['lambda_R']
            forward=-logmeanexp(-beta*delta*a)/beta;reverse=logmeanexp(beta*delta*b)/beta
            bins=np.linspace(min(a.min(),b.min()),max(a.max(),b.max())+1e-8,21)
            pa=np.histogram(a,bins)[0]/len(a);pb=np.histogram(b,bins)[0]/len(b)
            edge.update({'gap_histogram_overlap':float(np.minimum(pa,pb).sum()),'forward_deltaG_kj_mol':forward,'reverse_deltaG_kj_mol':reverse,'hysteresis_kj_mol':abs(forward-reverse)})
        edges.append(edge)
    result={'production_complete':raw['cursor']['phase']=='complete','sample_counts':[len(w['samples']) for w in windows],
        'edges':edges,'calibrated':False,'convergence_established':False,
        'note':'Histogram overlap and forward/reverse FEP are necessary checks, not proof of convergence. Samples may be correlated; no independent-replica uncertainty is inferred.'}
    atomic_json(out/'diagnostics.json',result);print(json.dumps(result,indent=2));return result


def run(args):
    cfg=load_config(args.config)
    if args.case is not None:cfg['case']=args.case
    out=args.output_dir or ROOT/'runs/full_enzyme_short';out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if args.mode=='diagnose':return diagnose(out,cfg['seed'],getattr(args,'records',None))
    if args.mode=='scan':
        from workshop_scan import run as run_scan
        return run_scan(args,out,cfg['seed'])
    if args.mode=='target-fit':
        from kemp_profile import run_target
        return run_target(args,out,cfg['seed'])
    if args.mode=='profile':
        from kemp_profile import run as run_profile
        return run_profile(args,out,cfg['seed'])
    signature,provenance=fingerprint(cfg)
    store=Store(out/'run.sqlite',signature,args.resume)
    evb,data,meta=build(cfg['case'],out)
    meta.update({'fingerprint':signature,'openmm_version':mm.version.version,'config_file':str(args.config) if args.config else None})
    atomic_json(out/'provenance.json',provenance);atomic_json(out/'model.json',meta)
    platform=mm.Platform.getPlatformByName(cfg['platform']);properties={'Threads':str(cfg['cpu_threads'])} if cfg['platform']=='CPU' else {}
    integrator=mm.LangevinMiddleIntegrator(cfg['temperature_kelvin'],cfg['friction_per_ps'],cfg['timestep_fs']*.001);integrator.setRandomNumberSeed(cfg['seed']);integrator.setConstraintTolerance(1e-6)
    sim=app.Simulation(evb.states['R'].topology,evb.evb_scaffold,integrator,platform,properties)
    sim.context.setPositions(np.array(data['coordinates_angstrom'])*.1);evb.set_state_lambdas(sim,[1.])
    store.restore(sim);cursor=store.cursor()
    started=time.monotonic();limit=args.max_runtime_seconds
    def expired():return limit is not None and time.monotonic()-started>=limit
    segment_start=sim.currentStep
    segment=f'{segment_start:012d}_{uuid.uuid4().hex[:10]}'
    sim.reporters.append(app.DCDReporter(str(out/f'trajectory_{segment}.dcd'),cfg['trajectory_interval']))
    trace_path=out/f'trace_{segment}.csv'
    sim.reporters.append(app.StateDataReporter(str(trace_path),cfg['checkpoint_interval'],step=True,time=True,potentialEnergy=True,kineticEnergy=True,temperature=True,speed=True))
    if not (out/'topology.pdb').exists():
        with (out/'topology.pdb').open('w') as handle:app.PDBFile.writeFile(sim.topology,sim.context.getState(getPositions=True).getPositions(),handle)
    def save(cursor,rows=(),export=False):
        state,_,positions,_=finite_state(sim)
        mobile=sum(sim.system.getParticleMass(i).value_in_unit(unit.dalton)>0 for i in range(sim.system.getNumParticles()))
        dof=3*mobile-sim.system.getNumConstraints()
        temperature=2*state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(dof*.00831446261815324)
        if temperature>cfg['max_temperature_kelvin']:raise ValueError(f'Instantaneous temperature {temperature:.1f} K exceeds the stability limit; refusing to checkpoint an unstable state.')
        if np.linalg.norm(np.ptp(positions,axis=0))>=9.9:raise ValueError('Droplet extent exceeds the direct-interaction cutoff; inspect escaping atoms.')
        cursor={**cursor,'step':sim.currentStep}
        store.save(sim,cursor,rows);atomic_json(out/'progress.json',cursor)
        if export:store.export(out,cfg,meta)
        print(f'{cursor["phase"]}: step {sim.currentStep}'+(f', window {min(cursor.get("window",0)+1,len(cfg['lambda_schedule']))}/{len(cfg["lambda_schedule"])}' if cursor['phase']=='production' else ''),flush=True)
    try:
        mode=args.mode
        if mode=='prepare':
            finite_state(sim)
            if cursor['phase']=='prepared':save(cursor)
            return
        if mode in ('run','relax') and cursor['phase'] in ('prepared','relaxing'):
            evb.set_state_lambdas(sim,[1.]);result=bounded_minimize(sim,cfg,out,expired)
            if result['stop_reason']=='time_limit':save({'phase':'relaxing'});return
            sim.context.setVelocitiesToTemperature(cfg['heating_temperatures_kelvin'][0],cfg['seed'])
            cursor={'phase':'relaxed'};save(cursor)
        if mode=='relax':return
        if mode in ('run','equilibrate') and cursor['phase'] in ('relaxed','equilibrating'):
            done=cursor.get('completed_steps',0);stages=[(t,cfg['heating_steps_per_temperature']) for t in cfg['heating_temperatures_kelvin']]+[(cfg['temperature_kelvin'],cfg['equilibration_steps'])]
            total=0
            for temperature,steps in stages:
                target=total+steps;integrator.setTemperature(temperature)
                while done<target:
                    if expired():return
                    n=min(cfg['checkpoint_interval'],target-done);sim.step(n);done+=n
                    cursor={'phase':'equilibrating','completed_steps':done,'temperature_kelvin':temperature};save(cursor)
                total=target
            cursor={'phase':'equilibrated'};save(cursor)
        if mode=='equilibrate':
            if cursor['phase'] not in ('equilibrated','production','complete'):raise ValueError('Run relaxation before equilibration.')
            return
        if mode not in ('run','production'):return
        if cursor['phase'] not in ('equilibrated','production','complete'):raise ValueError('Run relaxation and equilibration before production, or use --mode run.')
        if cursor['phase']=='complete':print('Production already complete.',flush=True);return
        integrator.setTemperature(cfg['temperature_kelvin'])
        for window in range(cursor.get('window',0),len(cfg['lambda_schedule'])):
            lam=cfg['lambda_schedule'][window];evb.set_state_lambdas(sim,[lam])
            equil_done=cursor.get('window_equilibration_done',0) if cursor.get('window')==window else 0
            prod_done=cursor.get('production_done',0) if cursor.get('window')==window else 0
            while equil_done<cfg['window_equilibration_steps']:
                if expired():return
                n=min(cfg['checkpoint_interval'],cfg['window_equilibration_steps']-equil_done);sim.step(n);equil_done+=n
                cursor={'phase':'production','window':window,'window_equilibration_done':equil_done,'production_done':prod_done};save(cursor)
            pending=[]
            while prod_done<cfg['production_steps_per_window']:
                if expired():
                    if pending:save(cursor,pending)
                    return
                n=min(cfg['sample_interval']-prod_done%cfg['sample_interval'],cfg['production_steps_per_window']-prod_done)
                sim.step(n);prod_done+=n
                rows=[]
                if prod_done%cfg['sample_interval']==0:
                    energies=tuple(float(v) for v in evb._sample_diabatic_state_energies_kj_mol(sim))
                    _,mapping,_,_=finite_state(sim)
                    if not all(math.isfinite(v) for v in energies):raise ValueError('Nonfinite sampled diabat energy.')
                    closure=mapping-(lam*energies[0]+(1-lam)*energies[1])
                    if abs(closure)>.05:raise ValueError(f'Diabat mapping closure failed: {closure} kJ/mol.')
                    rows=[(window,prod_done,sim.currentStep,lam,*energies,mapping)]
                cursor={'phase':'production','window':window,'window_equilibration_done':equil_done,'production_done':prod_done}
                pending.extend(rows)
                if prod_done%cfg['checkpoint_interval']==0 or prod_done==cfg['production_steps_per_window']:
                    save(cursor,pending);pending=[]
            cursor={'phase':'production','window':window+1,'window_equilibration_done':0,'production_done':0};save(cursor,export=True)
        cursor={'phase':'complete','windows_completed':len(cfg['lambda_schedule'])};save(cursor)
        print(f'Production finished: {out}. Raw mapping data are uncalibrated; no physical barrier is inferred.',flush=True)
    except (KeyboardInterrupt,Exception):
        # Last committed SQLite checkpoint and samples remain mutually consistent.
        store.export(out,cfg,meta)
        raise
    finally:
        store.export(out,cfg,meta)
        atomic_json(out/f'segment_{segment}.json',{'trajectory':f'trajectory_{segment}.dcd','trace':trace_path.name,'start_step':segment_start,'last_committed_step':store.cursor().get('step',segment_start),'context_step_at_exit':int(sim.context.getStepCount()),'note':'Frames beyond last_committed_step can overlap a resumed segment after an interruption.'})
        store.db.close()
