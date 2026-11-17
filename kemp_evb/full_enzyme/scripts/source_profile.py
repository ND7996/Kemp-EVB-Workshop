"""Archived QFEP-parameter comparison using per-frame EVB energy mixing.

This is not a new fit or a declaration of calibration for the changed OpenMM
boundary. Dynamics records are read without modifying the simulation or inputs.
"""
import hashlib
import json
import math
from pathlib import Path
import webbrowser
import numpy as np
from evbopenmm.models import FEPSampleRecord,FEPWindowRecord
from evbopenmm.analysis.mapping_fep import analyze_mapping_fep_umbrella_records
from evbopenmm.analysis.profile import detect_profile_thermochemistry
from kemp_inputs import ROOT

SOURCE=ROOT/'inputs/protocol_reference/qfep.inp'


def read_qfep(path=SOURCE):
    """Read the supplied two-state, one constant off-diagonal QFEP input strictly."""
    lines=[line.split('#',1)[0].strip() for line in Path(path).read_text().splitlines()]
    rows=[line.split() for line in lines if line]
    nfiles=int(rows[0][0]);states,predefined=map(int,rows[1][:2])
    if (states,predefined)!=(2,0):raise ValueError('Expected two states and externally specified coupling.')
    kt=float(rows[2][0]);skip=int(rows[2][1]);bins=int(rows[3][0]);minimum=int(rows[4][0]);alpha=float(rows[5][0])
    if int(rows[6][0])!=1:raise ValueError('Expected one coupling entry.')
    i,j=int(rows[7][0]),int(rows[7][1]);amplitude,mu,eta,r0=map(float,rows[7][2:6])
    if (i,j)!=(1,2) or any(v!=0 for v in (mu,eta,r0)):raise ValueError('Archived coupling is not the supported constant form.')
    if list(map(float,rows[8]))!=[1.,-1.]:raise ValueError('Expected reaction coordinate state1-state2.')
    files=[row[0] for row in rows[9:] if row[0]!='stop']
    if len(files)!=nfiles or any(not name.endswith('.en') for name in files):raise ValueError('QFEP file list differs from declared count.')
    if not all(math.isfinite(v) for v in (kt,alpha,amplitude)) or kt<=0 or amplitude<=0:raise ValueError('Invalid QFEP energy parameters.')
    return {'source_dataset_doi':'10.5281/zenodo.12666978','qfep_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        'state_order':['R','P'],'native_alpha_kcal_mol':alpha,'coupling_kcal_mol':amplitude,'coupling_form':'constant',
        'native_kT_kcal_mol':kt,'native_skip_samples_per_window':skip,'native_bins':bins,'native_min_points_per_window_bin':minimum,
        'native_windows':nfiles,'energy_filenames':files,'transfer_validated_for_current_model':False}


def records(raw,skip=0):
    if raw.get('schema')!='kemp_full_mapping_records_v1':raise ValueError('Expected full-enzyme raw mapping records.')
    if skip<0:raise ValueError('skip_samples must be nonnegative.')
    result=[]
    for window in raw['windows']:
        lam=float(window['lambda_R']);samples=[]
        if not math.isfinite(lam) or not 0<=lam<=1:raise ValueError('Invalid mapping weight.')
        for row in window['samples'][skip:]:
            r,p=float(row['E_R_kj_mol']),float(row['E_P_kj_mol'])
            if not all(math.isfinite(v) for v in (r,p)):raise ValueError('Nonfinite sampled energy.')
            samples.append(FEPSampleRecord(int(row['step']),lam*r+(1-lam)*p,r-p,(r,p)))
        if not samples:raise ValueError('No samples remain in a window; native skip4000 cannot be applied to this short run.')
        result.append(FEPWindowRecord(int(window['index']),(lam,1-lam),tuple(samples)))
    if len(result)<2:raise ValueError('At least two sampled mapping windows are required.')
    return tuple(result)


def frame_bin_free_energy(deltas,kT,window_free_energy):
    """Conditional-bin exponential average of per-frame target-minus-mapping E."""
    values=np.asarray(deltas,dtype=float)
    if not len(values) or not np.isfinite(values).all() or kT<=0:raise ValueError('Finite nonempty bin data and positive kT required.')
    exponents=-values/kT;maximum=float(exponents.max())
    return float(window_free_energy-kT*(maximum+np.log(np.exp(exponents-maximum).mean())))


def reconstruct_source(raw,bins=40,min_points=20,skip_samples=0,parameters=None):
    if bins<3 or min_points<1:raise ValueError('Use at least three bins and one point per window/bin.')
    params=read_qfep() if parameters is None else parameters
    windows=records(raw,skip_samples)
    temperature=float(raw['config']['temperature_kelvin'])
    if not math.isfinite(temperature) or temperature<=0:raise ValueError('Invalid sampling temperature.')
    kt=.00831446261815324*temperature/4.184
    correction=float(raw['metadata']['morse_P_minus_R_correction_kcal_mol'])
    alpha=float(params['native_alpha_kcal_mol']);offset=alpha+correction;coupling=float(params['coupling_kcal_mol'])
    if not all(math.isfinite(v) for v in (correction,alpha,coupling)) or coupling<=0:raise ValueError('Invalid source comparison parameters.')
    base=analyze_mapping_fep_umbrella_records(windows,state_names=('R','P'),kT_kj_mol=kt*4.184,skip_samples=0,n_gap_bins=bins,min_points_per_bin=min_points,input_energy_unit='kJ/mol')
    gap=np.array([s.reaction_coordinate_kj_mol/4.184 for w in windows for s in w.samples])
    if gap.max()<=gap.min():raise ValueError('The sampled gap has no range.')
    edges=np.linspace(gap.min(),np.nextafter(gap.max(),np.inf),bins+1)
    weighted=np.zeros(bins);accepted=np.zeros(bins,dtype=int);raw_counts=np.histogram(gap,edges)[0]
    for index,window in enumerate(windows):
        g=np.array([s.reaction_coordinate_kj_mol/4.184 for s in window.samples]);lam=window.lambda_values[0]
        # Subtract the common P energy analytically for numerical stability.
        ground_relative_P=.5*(g+offset)-np.hypot(.5*(g-offset),coupling)
        target_minus_mapping=ground_relative_P-lam*g
        assignments=np.clip(np.searchsorted(edges,g,side='right')-1,0,bins-1)
        for b in range(bins):
            values=target_minus_mapping[assignments==b];n=len(values)
            if n<min_points:continue
            estimate=frame_bin_free_energy(values,kt,base['fep']['cumulative_average_kcal_mol'][index])
            weighted[b]+=n*estimate;accepted[b]+=n
    curve=[float(weighted[b]/accepted[b]) if accepted[b] else None for b in range(bins)]
    centers=(edges[:-1]+edges[1:])*.5-offset  # QFEP uses shifted E_R-E_P.
    valid=[i for i,y in enumerate(curve) if y is not None]
    missing=[i for i in range(min(valid),max(valid)+1) if i not in valid] if valid else []
    thermo={'reactant':None,'product':None,'transition_states':(),'warnings':[]}
    if len(valid)>=3 and not missing:
        local=detect_profile_thermochemistry([float(centers[i]) for i in valid],[curve[i] for i in valid]);thermo.update(local)
        for key in ('reactant','product'):
            if thermo[key] is not None:thermo[key]={**thermo[key],'index':valid[int(thermo[key]['index'])]}
        thermo['transition_states']=tuple({**item,'index':valid[int(item['index'])]} for item in local['transition_states'])
    barrier=None;dg=None
    if thermo['reactant'] is not None and thermo['product'] is not None:
        dg=thermo['product']['free_energy_kcal_mol']-thermo['reactant']['free_energy_kcal_mol']
        if thermo['transition_states']:barrier=max(item['free_energy_kcal_mol'] for item in thermo['transition_states'])-thermo['reactant']['free_energy_kcal_mol']
    warnings=list(thermo['warnings'])
    warnings.extend(['Archived GNCA QFEP parameters, not newly fitted to13.2 or validated for this changed OpenMM boundary.',
        'Short correlated single-replica samples do not establish convergence; no confidence interval is inferred.',
        'Per-frame mixing before exponential averaging differs from the earlier illustrative eigenvalue-of-binned-free-energies construction.',
        'Native binning, burn-in and minimum counts are recorded but not copied blindly into the short-run comparison.'])
    if len(valid)<3:warnings.append('Fewer than three bins meet the per-window minimum; no barrier can be identified.')
    if missing:warnings.append('Internal unresolved bins: no barrier is reported across the disconnected profile.')
    if skip_samples!=params['native_skip_samples_per_window']:warnings.append(f'Actual extra skip={skip_samples}; archived analysis uses {params["native_skip_samples_per_window"]} records/window. Original simulation equilibration was discarded before these raw production records.')
    analysis={'profile':{'gap_centers_kcal_mol':tuple(float(x) for x in centers),'state_free_energies_kcal_mol':{'evb_ground':tuple(curve)},'bin_counts':tuple(int(v) for v in accepted)},'thermochemistry':thermo}
    summary={'analysis_method':'QFEP-style per-frame EVB mixing and within-window/bin exponential reweighting; count-weighted window combination',
        'calibrated':False,'physical_barrier_available':False,'convergence_established':False,'author_parameters_recovered':True,
        'coupling_kcal_mol':coupling,'extra_product_offset_kcal_mol':alpha,'morse_convention_correction_kcal_mol':correction,'total_product_offset_kcal_mol':offset,
        'conditional_barrier_estimate_kcal_mol':barrier,'conditional_reaction_free_energy_kcal_mol':dg,'resolved_bins':len(valid),'internal_missing_bins':missing,
        'reference':'reactant' if thermo['reactant'] is not None else 'minimum','plot_title':'Kemp EVB — ARCHIVED PARAMETERS; TRANSFER NOT VALIDATED',
        'gap_axis_label':'Energy gap (kcal/mol)','gap_csv_column':'offset_adjusted_energy_gap_kcal_mol',
        'sampling_complete':raw['cursor']['phase']=='complete','sample_count':sum(len(w.samples) for w in windows),'sampled_windows':len(windows),
        'bins':bins,'min_points_per_window_bin':min_points,'actual_extra_skip_samples_per_window':skip_samples,'raw_bin_counts':[int(v) for v in raw_counts],
        'accepted_bin_counts':[int(v) for v in accepted],'temperature_kelvin':temperature,'native_qfep':params,'thermochemistry':thermo,
        'not_identical_to_native_QFEP':'Different raw energies/boundary/sampling, configurable binning and skip; highest-gap endpoint is retained explicitly. No native binary .en parity run has been performed.',
        'warnings':warnings}
    summary['case']=raw.get('metadata',{}).get('case') or raw.get('config',{}).get('case')
    summary['parameter_provenance_note']='Coupling/offset recovered from the separate GNCA4-WT QFEP reference archive (Zenodo 12666978), not verified original Des27.7-IN Q6 settings. Transfer to this OpenMM adaptation is unvalidated.'
    return analysis,summary


def run(args,out,seed):
    from kemp_profile import write_outputs
    path=getattr(args,'records',None) or Path(out)/f'replica_{seed}_raw.json';raw=json.loads(path.read_text())
    archived=read_qfep();params=archived
    if args.coupling is not None:
        params={**archived,'coupling_kcal_mol':args.coupling,'native_alpha_kcal_mol':args.offset}
    analysis,summary=reconstruct_source(raw,args.bins,args.min_points,args.skip_samples,params)
    summary['raw_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    if args.coupling is None:
        name=f'profile_source_qfep_bins{args.bins}_min{args.min_points}_skip{args.skip_samples}'
    else:
        name=f'profile_trial_H{args.coupling:g}_offset{args.offset:g}_bins{args.bins}_min{args.min_points}'
        summary.update({'author_parameters_recovered':False,'native_qfep':archived,'parameter_origin':'User-selected trial; not physical calibration.',
                        'plot_title':'Kemp EVB — USER PARAMETER TRIAL; NOT A PREDICTION'})
        summary['warnings'].insert(0,'User-selected coupling/offset trial, not the archived parameter set.')
    target=Path(out)/name;write_outputs(analysis,summary,target)
    print(f'Full-enzyme parameter-conditional profile: {target/"profile.html"}')
    print(f'H12={summary["coupling_kcal_mol"]:g}; extra product offset={summary["extra_product_offset_kcal_mol"]:g}; Morse correction={summary["morse_convention_correction_kcal_mol"]:+g} kcal/mol.')
    print(f'Conditional forward barrier: {summary["conditional_barrier_estimate_kcal_mol"]} kcal/mol. Not a validated physical barrier.')
    if summary['internal_missing_bins']:print('Warning: unresolved internal bins; no barrier is inferred across a gap.')
    if args.open:webbrowser.open((target/'profile.html').resolve().as_uri())
